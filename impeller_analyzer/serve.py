"""Application locale : deposer une roue dans le navigateur, l'analyser sur sa machine.

    python -m impeller_analyzer.serve

Le calcul est en Python et ne peut pas s'executer dans le navigateur.  Plutot
que d'en tenir une seconde version en JavaScript -- deux implementations de la
meme physique a garder d'accord, dont une sans tests -- ce module sert la page
et fait tourner **le meme code que la ligne de commande** derriere un petit
serveur HTTP de la bibliotheque standard.

Le serveur n'ecoute que sur la machine locale et ne parle a personne d'autre :
le fichier depose ne quitte pas le poste.
"""

from __future__ import annotations

import argparse
import http.server
import json
import os
import shutil
import socketserver
import sys
import tempfile
import threading
import traceback
import urllib.parse
import uuid
import webbrowser

from . import config
from .analysis import Options, run
from .geometry import blade_angles
from .io import loader, report, viewer

#: Extensions acceptees par le formulaire, toutes familles confondues.
ACCEPTED = (
    config.EXT_MESH_NATIVE
    + config.EXT_MESH_TRIMESH_ONLY
    + config.EXT_CAD_BREP
    + config.EXT_DXF
    + config.EXT_DWG
)

#: Fichiers qu'une analyse met a disposition, dans l'ordre d'affichage.
DOWNLOADABLE = (
    report.MARKDOWN_NAME,
    report.JSON_NAME,
    report.CURVES_NAME,
    report.OCCUPANCY_NAME,
    report.VIEWER_NAME,
)


class Store:
    """Repertoires d'analyse, gardes le temps de telecharger les rapports."""

    def __init__(self, root: str, keep: int = config.SERVER_HISTORY):
        self.root = root
        self.keep = keep
        self.order: list[str] = []
        self.paths: dict[str, str] = {}
        self.lock = threading.Lock()

    def create(self) -> tuple[str, str]:
        """Nouveau repertoire d'analyse ; renvoie `(identifiant, chemin)`."""
        token = uuid.uuid4().hex
        path = os.path.join(self.root, token)
        os.makedirs(path, exist_ok=True)
        with self.lock:
            self.paths[token] = path
            self.order.append(token)
            while len(self.order) > self.keep:
                old = self.order.pop(0)
                shutil.rmtree(self.paths.pop(old, ""), ignore_errors=True)
        return token, path

    def resolve(self, token: str, name: str) -> str | None:
        """Chemin d'un fichier produit, ou `None` si la demande n'est pas legitime."""
        if name not in DOWNLOADABLE:
            return None
        with self.lock:
            base = self.paths.get(token)
        if base is None:
            return None
        candidate = os.path.join(base, name)
        return candidate if os.path.isfile(candidate) else None

    def clear(self) -> None:
        """Efface tous les repertoires d'analyse."""
        with self.lock:
            for path in self.paths.values():
                shutil.rmtree(path, ignore_errors=True)
            self.paths.clear()
            self.order.clear()


class BadRequest(Exception):
    """Demande invalide, a renvoyer telle quelle a l'utilisateur."""


def _float(query: dict, key: str, default=None):
    """Lit un nombre du formulaire, avec un message clair s'il est illisible."""
    raw = query.get(key, [""])[0].strip().replace(",", ".")
    if not raw:
        return default
    try:
        return float(raw)
    except ValueError:
        raise BadRequest(f"« {key} » doit etre un nombre, or « {raw} » n'en est pas un") from None


def _int(query: dict, key: str, default=None):
    """Lit un entier du formulaire."""
    value = _float(query, key)
    if value is None:
        return default
    return int(round(value))


def parse_speeds(raw: str) -> tuple[float, ...]:
    """Lit une liste de regimes separes par des espaces ou des virgules."""
    tokens = [token for token in raw.replace(",", " ").split() if token]
    if not tokens:
        return config.DEFAULT_RPM
    speeds = []
    for token in tokens:
        try:
            value = float(token)
        except ValueError:
            raise BadRequest(f"regime illisible : « {token} »") from None
        if value <= 0.0:
            raise BadRequest("les regimes doivent etre strictement positifs")
        speeds.append(value)
    return tuple(speeds)


def _rotation(query: dict) -> int | None:
    """Sens de rotation demande par le formulaire, `None` s'il n'est pas renseigne."""
    value = (query.get("rotation", [""])[0] or "").strip()
    if not value:
        return None
    try:
        return blade_angles.rotation_sign_from_name(value)
    except ValueError as error:
        raise BadRequest(str(error)) from error


def options_from_query(query: dict) -> Options:
    """Construit les options d'analyse depuis les champs du formulaire."""
    grid = _int(query, "grille", config.GRID_NR) or config.GRID_NR
    grid = max(40, min(400, grid))
    return Options(
        unit=query.get("unite", ["cm"])[0] or "cm",
        r_aspiration_cm=_float(query, "r_aspiration"),
        speeds=parse_speeds(query.get("rpm", [""])[0]),
        blades=_int(query, "pales"),
        beta1_deg=_float(query, "beta1"),
        beta2_deg=_float(query, "beta2"),
        rotation=_rotation(query),
        altitude=_float(query, "altitude", config.ALTITUDE),
        temperature_c=_float(query, "temperature", config.TEMPERATURE),
        suction_height=_float(query, "hauteur", config.HAUTEUR_ASPIRATION),
        suction_losses=_float(query, "pertes", config.PERTES_ASPIRATION),
        grid_nr=grid,
        grid_nz=grid,
    )


def extension_of(name: str) -> str:
    """Extension validee du fichier depose."""
    extension = os.path.splitext(os.path.basename(name))[1].lower()
    if extension in config.EXT_REFUSED:
        raise BadRequest(
            f"un fichier {extension} est un programme AutoLISP, pas un format geometrique. "
            "Exportez la roue en STL depuis votre logiciel de CAO."
        )
    if extension not in ACCEPTED:
        raise BadRequest(
            f"extension non geree : {extension or '(aucune)'}. Formats lus : {', '.join(ACCEPTED)}"
        )
    return extension


def analyse_upload(payload: bytes, query: dict, store: Store) -> dict:
    """Analyse un fichier depose et renvoie la charge utile de la vue 3D."""
    name = os.path.basename(query.get("nom", ["roue.stl"])[0]) or "roue.stl"
    extension = extension_of(name)
    if not payload:
        raise BadRequest("fichier vide")

    options = options_from_query(query)
    token, directory = store.create()
    source = os.path.join(directory, "entree" + extension)
    with open(source, "wb") as handle:
        handle.write(payload)

    try:
        result = run(source, options)
    except loader.ImportError_ as error:
        raise BadRequest(str(error)) from None
    except ValueError as error:
        raise BadRequest(str(error)) from None

    if result.import_report is not None:
        result.import_report.path = name  # le nom d'origine, pas le fichier temporaire
    produced = report.write_all(result, directory, source=name)
    data = viewer.build_payload(result)
    data["telechargements"] = [
        {"nom": os.path.basename(path), "url": f"/telecharger/{token}/{os.path.basename(path)}"}
        for key, path in produced.items()
        if os.path.basename(path) in DOWNLOADABLE
    ]
    data["telechargements"].sort(key=lambda item: DOWNLOADABLE.index(item["nom"]))
    return data


class Handler(http.server.BaseHTTPRequestHandler):
    """Trois routes : la page, l'analyse, le telechargement des rapports."""

    server_version = "impeller-analyzer"
    store: Store
    page: str

    def log_message(self, fmt, *args):  # pragma: no cover - bruit de console
        """Journal reduit : une ligne par requete, sans l'horodatage complet."""
        sys.stderr.write("  %s\n" % (fmt % args))

    # -- envois -----------------------------------------------------------
    def _send(self, code: int, body: bytes, content_type: str, extra: dict | None = None) -> None:
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for key, value in (extra or {}).items():
            self.send_header(key, value)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, code: int, data: dict) -> None:
        payload = json.dumps(data, ensure_ascii=False, allow_nan=False).encode("utf-8")
        self._send(code, payload, "application/json; charset=utf-8")

    # -- routes -----------------------------------------------------------
    def do_GET(self) -> None:  # noqa: N802 - nom impose par BaseHTTPRequestHandler
        """Sert la page, un fichier produit, ou 404."""
        route = urllib.parse.urlparse(self.path).path
        if route in ("/", "/index.html"):
            self._send(200, self.page.encode("utf-8"), "text/html; charset=utf-8")
            return
        parts = [part for part in route.split("/") if part]
        if len(parts) == 3 and parts[0] == "telecharger":
            path = self.store.resolve(parts[1], parts[2])
            if path is None:
                self._send(404, b"fichier introuvable", "text/plain; charset=utf-8")
                return
            with open(path, "rb") as handle:
                body = handle.read()
            kind = "text/html; charset=utf-8" if path.endswith(".html") else (
                "image/png" if path.endswith(".png") else (
                    "application/json; charset=utf-8" if path.endswith(".json")
                    else "text/markdown; charset=utf-8"))
            self._send(200, body, kind,
                       {"Content-Disposition": f'attachment; filename="{parts[2]}"'})
            return
        self._send(404, b"page introuvable", "text/plain; charset=utf-8")

    def do_HEAD(self) -> None:  # noqa: N802
        """Meme routage que GET, sans corps."""
        self.do_GET()

    def do_POST(self) -> None:  # noqa: N802
        """Recoit le fichier depose et renvoie le resultat de l'analyse."""
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path != "/analyse":
            self._json(404, {"erreur": "route inconnue"})
            return
        try:
            length = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            self._json(400, {"erreur": "en-tete Content-Length illisible"})
            return
        if length <= 0:
            self._json(400, {"erreur": "aucun fichier recu"})
            return
        if length > config.SERVER_MAX_UPLOAD:
            self._json(413, {
                "erreur": f"fichier trop volumineux ({length / 1048576:.0f} Mo, "
                          f"maximum {config.SERVER_MAX_UPLOAD / 1048576:.0f} Mo)"
            })
            return

        body = self.rfile.read(length)
        query = urllib.parse.parse_qs(parsed.query, keep_blank_values=True)
        try:
            self._json(200, analyse_upload(body, query, self.store))
        except BadRequest as error:
            self._json(400, {"erreur": str(error)})
        except Exception:  # pragma: no cover - filet de securite
            traceback.print_exc()
            self._json(500, {"erreur": "l'analyse a echoue ; la trace est dans la console du serveur"})


class Server(socketserver.ThreadingMixIn, http.server.HTTPServer):
    """Serveur local, une requete a la fois par connexion."""

    daemon_threads = True
    allow_reuse_address = True


def build_server(host: str = config.SERVER_HOST, port: int = config.SERVER_PORT) -> tuple[Server, Store]:
    """Prepare le serveur et son magasin d'analyses."""
    store = Store(tempfile.mkdtemp(prefix="impeller-app-"))
    handler = type("BoundHandler", (Handler,), {"store": store, "page": viewer.build_app_page()})
    return Server((host, port), handler), store


def main(argv: list[str] | None = None) -> int:  # pragma: no cover - point d'entree
    """Lance l'application locale."""
    parser = argparse.ArgumentParser(
        prog="python -m impeller_analyzer.serve",
        description=(
            "Application locale : deposez une roue dans le navigateur, l'analyse s'execute "
            "sur votre machine par le meme code que la ligne de commande."
        ),
    )
    parser.add_argument("--port", type=int, default=config.SERVER_PORT, help="port d'ecoute (defaut : %(default)s)")
    parser.add_argument(
        "--host", default=config.SERVER_HOST,
        help="adresse d'ecoute (defaut : %(default)s, c'est-a-dire cette machine seulement)",
    )
    parser.add_argument("--sans-navigateur", action="store_true", help="ne pas ouvrir le navigateur")
    args = parser.parse_args(argv)

    try:
        server, store = build_server(args.host, args.port)
    except OSError as error:
        print(f"impossible d'ecouter sur {args.host}:{args.port} : {error}", file=sys.stderr)
        print("le port est peut-etre deja pris ; essayez --port 8766", file=sys.stderr)
        return 2

    url = f"http://{args.host}:{args.port}/"
    print(f"Inspecteur de roue : {url}")
    print("  deposez un fichier 3D dans la page ; le calcul se fait ici, rien ne sort de la machine")
    if args.host not in ("127.0.0.1", "localhost", "::1"):
        print(f"  attention : l'ecoute sur {args.host} rend l'application accessible depuis le reseau")
    print("  Ctrl+C pour arreter")
    if not args.sans_navigateur:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\narret")
    finally:
        server.server_close()
        store.clear()
        shutil.rmtree(store.root, ignore_errors=True)
    return 0


if __name__ == "__main__":  # pragma: no cover - point d'entree
    sys.exit(main())
