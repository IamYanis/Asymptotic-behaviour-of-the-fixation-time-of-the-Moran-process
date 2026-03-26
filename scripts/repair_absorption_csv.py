# repair_absorption_csv.py
from pathlib import Path
import re
import sys
import pandas as pd

src = Path("data/absorption_times.csv")
dst = Path("data/absorption_times.cleaned.csv")

if not src.exists():
    sys.exit(f"Fichier introuvable: {src}")

# Remplace uniquement la virgule entre deux chiffres (ex: 0,123 -> 0.123)
decimal_comma = re.compile(r'(?<=\d),(?=\d)')

lines_fixed = 0
with src.open("r", encoding="utf-8", newline="") as fin, \
     dst.open("w", encoding="utf-8", newline="") as fout:
    header = fin.readline()
    fout.write(header)
    expected_cols = len([c.strip() for c in header.rstrip("\n\r").split(",")])
    for i, line in enumerate(fin, start=2):
        new_line = decimal_comma.sub(".", line)
        if new_line != line:
            lines_fixed += 1
        # Optionnel: on peut enlever les espaces autour des virgules
        # new_line = re.sub(r'\s*,\s*', ',', new_line)
        fout.write(new_line)

print(f"Remplacements virgule->point effectués sur {lines_fixed} lignes.")
print(f"Fichier nettoyé écrit: {dst}")

# Vérifier la lecture et le nombre de colonnes (attendu: 7)
try:
    df = pd.read_csv(dst)
    print(f"Lecture OK. Colonnes: {list(df.columns)} (n={len(df)})")
    if len(df.columns) != 7:
        print("⚠️ Avertissement: nombre de colonnes inattendu (attendu 7).")
except Exception as e:
    print("Échec de lecture du fichier nettoyé avec pandas.")
    print(e)
    print("Astuce: si le fichier a été ré-enregistré par Excel en ';' et décimales ',', essaye:")
    print("pd.read_csv('data/absorption_times.csv', sep=';', decimal=',')")
