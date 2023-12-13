<h2 align="center">phylogenetic tree - UPGMA</h2>

<p align="center">Árbol filogenético con el método de distancia UPGMA</p>

## 📍 Descripción

Este proyecto es un script en línea de comandos para generar un árbol filogenético empleando el método de distancia "Unweighted Pair Group Method with Arithmetic mean" (UPGMA). Este es un método de clustering que se basa en la identificación de las parejas más similares y en el cálculo de la media de las distancias entre ellas y el resto de las secuencias para reconstruir el árbol.

Emplea este script para generar un árbol.

## 📌 Prerequisitos

- Python 3.0 > [instalar](https://www.python.org/)
- biopython 1.81 > [instalar](https://biopython.org/wiki/Download)
- ete3 3.1.3 > [instalar](https://pypi.org/project/ete3/)

## 📌 Instalación

```
git clone https://github.com/jmluna123/bioScript.git

cd .\UPGMA

pip install -r requirements.txt
```

## 📌 Empezar

El script solo requiere la ruta del archivo .fasta que contenga las diferentes secuencias a comparar.

El script muestra el árbol generado usando el siguiente comando:

```
python .\UPGMA.py <ruta archivo FASTA>
```

## 📌 Ejemplos

Para obtener y visualizar el árbol:

```
python .\UPGMA.py .\Crab_rRNA.fasta
```
