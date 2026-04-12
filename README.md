# Risolutore-FEM-python
Fem-travi è un risolutore bidimensionale basato sul metedo degli elementi finiti (FEM).

Il progetto nasce dall'esigenza pratica di avere uno strumento di visualizzazione per la risoluzione di esercizi e prove d'esame di meccanica dei solidi. Lo strumento permette di verificare rapidamente spostamenti, reazioni vincolari e diagrammi delle sollecitazioni iperstatiche.


# Funzionalità principali
Modellazione fisica: trave eulero bernulli con 3 gdl per nodo (ux, uy,\phi$) 
Carichi complessi: supporto per forze nodali, momenti e carichi distribuiti (uniformi e triangolari)
Analisi iperstatica: Risoluzione del sistema globale tramite matrici di rigidezza ruotate e assemblate
Verifica usata: Calcolo delle reazioni vincolari e successivamente verifica dell'equilibrio)
Interfaccia grafica (GUI): Strumento visivo interattivo per la costruizione del modello e la verifica




# Installazione
Il progetto è pensato come un pacchetto python installabile
1. clona la repository
   git clone https://github.com/Mariglend/Risolutore-FEM-python.git
   cd fem-travi

2. installa pacchetto
   pip install -e ".[dev]"




#Utilizzo
Il progetto è pensato per essere utilizzato sia come libreria python che tramite l'interfaccia dedicata.
Per avviare l'interfaccia da terminale eseguire:
python gui_fem.py





#Script rapido
Per verifiche veloci
from fem_travi import Struttura, Nodo, Trave

s = Struttura()
r = s.risolvi()
r.stampa()





# Validazione e Metodologia
Il solver utilizza il metodo della rigidezza diretta. La validazione è stata fatta confrontando i risultati con soluzioni analitiche di prove d'esame che coprono:

   
  1.Trave appoggiata a sbalzo
  
  2.Telai e portali iperstatici
  
  3.VErifica della simmetria e della deficinizia positiva della matrice K globale






OSS: Questo progetto è stato sviluppato per scopi didattici sebbene i risultati siano stati valididati si consiglia un controllo critico dei risultati.
