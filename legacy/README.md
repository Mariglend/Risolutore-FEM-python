# Versione precedente (v1)

- `gui_fem_pro_v2.py` — interfaccia tkinter della v1 (vincoli, aste, database sezioni/materiali, analisi modale, carichi termici).
  Usa un solver interno separato dalla libreria, con i due problemi corretti nella v2:
  il diagramma del momento sotto carichi distribuiti è lineare (dovrebbe essere parabolico)
  e le capriate di sole aste risultano labili.
- `solver_travi_fem.html` — prima versione web.
- `README_libreria_v1.md` — documentazione della libreria v1.

Restano qui come riferimento; la versione attuale è l'interfaccia `fem-travi`.
