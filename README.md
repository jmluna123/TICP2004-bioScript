<h2 align="center">BioScript</h2>

<p align="center">Procesos de la información genética</p>

## 📍 Descripción

Este proyecto es un script en línea de comandos para realizar procesos de transcripción y traducción de la información genética.
Estos procesos son los encargados de transformar la información en el ADN a proteínas, las cuales son las moléculas encargadas de realizar las diversas funcionalidades biológicas en nuestras células.

Emplea este script para transcribir y traducir sencuencias de forma rápida.

## 📌 Prerequisitos

- Python 3.0 > [instalar](https://www.python.org/)

## 📌 Instalación

```
git clone https://github.com/jmluna123/bioScript.git
```

## 📌 Empezar

El script cuenta con diferentes modos o flags para dar libertad al usuario en generar los archivos.

| Flag corta | Flag larga        | Por Defecto | Descripción                                                                                                                                           |
| ---------- | ----------------- | ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-h`       | `--help`          | -           | Muestra mensaje de ayuda con las opciones del script.                                                                                                 |
| `-v`       | `--verbose`       | `False`     | Muestra mensajes de los pasos que se efectuan al procesar las secuencias. Por defecto solo muestra mensajes de error, advertencia y archivos finales. |
| `-i`       | `--inverse`       | `False`     | Genera archivos de la secuencia inversa (3' a 5')                                                                                                     |
| `-tc`      | `--transcription` | `True`      | Genera archivos del proceso de transcripción o retrotranscripción de la secuencia ingresada (5' a 3')                                                 |
| `-tl`      | `--translation`   | `True`      | Genera archivos del proceso de traducción. Genera frame 1, 2 y 3.                                                                                     |

El script por defecto ejecuta hasta el proceso de traducción, usando el siguiente comando:

```
python .\bioScript.py <ruta archivo FASTA> [ -h help | -v verbose | -tc transcription | -rtc retrotranscription | -tl translation ]
```

## 📌 Ejemplos

Para obtener los archivos de **transcripción** (5' a 3') y **traducción** frame 1, 2 y 3:

```
python .\bioScript.py .\examples\example.fasta
```

Para obtener los archivos de los procesos y ver **mensajes del procesamiento**:

```
python .\bioScript.py .\examples\example.fasta -v
```

Para obtener solo los archivos de **transcripción** (5' a 3'):

```
python .\bioScript.py .\examples\tarea_transcribir.fasta -tc
```

Para obtener solo los archivos de **transcripción** (5' a 3') y (3' a 5'):

```
python .\bioScript.py .\examples\tarea_transcribir.fasta -tc -i
```

Para obtener solo los archivos de **traducción** frame 1, 2 y 3:

```
python .\bioScript.py .\examples\tarea_traducir.fasta -tl
```
