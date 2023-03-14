# bioinfo
Ricostruzione di una sequenza a partire dalle read attraverso un grafo di de Bruijn.

## Obiettivo
L'idea è di implementare l'algoritmo di ricostruzione di una sequenza a partire dalle read attraverso un grafo di de Bruijn. I passi sarebbero i seguenti:
- Scaricare un file contenente le read da https://www.ncbi.nlm.nih.gov/sra
- Processare il file con Biopython
- Costruire il grafo di de Bruijn
- Visualizzare il grafo con https://github.com/cytoscape/py4cytoscape 
- Ricostruire la sequenza
- Stampare la sequenza ricostruita

## Setup
Download, install and run Cytoscape: https://cytoscape.org/  
pip install python-igraph requests pandas networkx colorbrewer chardet decorator backoff colour  
pip install py4cytoscape  
pip install biopython  
