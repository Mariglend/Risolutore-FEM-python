"""
fem_travi/app/server.py
=======================
Interfaccia grafica: un piccolo server locale (solo libreria standard) che
serve la pagina web e fa i calcoli con la libreria fem_travi.

    fem-travi                 # apre il browser
    fem-travi --porta 9000 --no-browser
"""

from __future__ import annotations

import argparse
import json
import socket
import threading
import webbrowser
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from ..io import esporta_csv, risolvi_dict, struttura_da_dict

STATIC = Path(__file__).with_name("static")
MAX_BODY = 2_000_000


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=str(STATIC), **kw)

    def log_message(self, fmt, *args):  # silenzioso
        pass

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def _json_body(self):
        n = int(self.headers.get("Content-Length") or 0)
        if n <= 0 or n > MAX_BODY:
            raise ValueError("Richiesta vuota o troppo grande")
        return json.loads(self.rfile.read(n).decode("utf-8"))

    def _send(self, code, body: bytes, ctype: str):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        try:
            modello = self._json_body()
            if self.path == "/api/risolvi":
                out = risolvi_dict(modello, n_punti=int(modello.get("_punti", 41)))
                self._send(200, json.dumps(out).encode(), "application/json")
            elif self.path == "/api/csv":
                r = struttura_da_dict(modello).risolvi()
                self._send(200, esporta_csv(r).encode("utf-8-sig"), "text/csv; charset=utf-8")
            else:
                self._send(404, b"{}", "application/json")
        except Exception as e:  # l'interfaccia mostra il messaggio
            out = {"ok": False, "errore": f"{type(e).__name__}: {e}"}
            self._send(HTTPStatus.OK, json.dumps(out).encode(), "application/json")


def _porta_libera(preferita: int) -> int:
    for p in [preferita] + list(range(preferita + 1, preferita + 50)):
        with socket.socket() as s:
            try:
                s.bind(("127.0.0.1", p))
                return p
            except OSError:
                continue
    return 0


def avvia(porta: int = 8765, apri_browser: bool = True) -> None:
    porta = _porta_libera(porta)
    srv = ThreadingHTTPServer(("127.0.0.1", porta), Handler)
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    print(f"fem-travi: interfaccia su {url}  (Ctrl+C per chiudere)")
    if apri_browser:
        threading.Timer(0.4, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nChiuso.")
    finally:
        srv.server_close()


def main(argv=None):
    p = argparse.ArgumentParser(prog="fem-travi", description="Solver FEM travi 2D — interfaccia grafica")
    p.add_argument("--porta", type=int, default=8765)
    p.add_argument("--no-browser", action="store_true", help="non aprire il browser")
    a = p.parse_args(argv)
    avvia(a.porta, not a.no_browser)


if __name__ == "__main__":
    main()
