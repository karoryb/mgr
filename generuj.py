import os
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
 
def generuj_cgr(sekwencja):
    n = len(sekwencja)
    punkty = np.zeros((2, n+1))
    punkty[:, 0] = np.array([0.5, 0.5])  # start w środku

    for i in range(1, n+1):
        nukleotyd = sekwencja[i-1].upper()
        if nukleotyd == 'A':
            rog = np.array([0, 0])
        elif nukleotyd == 'C':
            rog = np.array([0, 1])
        elif nukleotyd == 'G':
            rog = np.array([1, 1])
        elif nukleotyd == 'T':
            rog = np.array([1, 0])
        else:
            continue  # gdzy jakis blad X nieznany itd

        punkty[:, i] = 0.5 * (punkty[:, i-1] + rog)

    return punkty

# # Wczytanie sekwencji z pliku FASTA
# plik_fasta = r"C:\Users\Karola\OneDrive - MULTIKKA Sp. z o.o\Dokumenty\STUDIA\INS\magisterka\NC_001700.1.fasta"

# for rekord in SeqIO.parse(plik_fasta, "fasta"):
#     sekwencja = str(rekord.seq)
#     punkty = generuj_cgr(sekwencja)
    
#     # opcjonalnie: rysowanie CGR
#     plt.figure(figsize=(5,5))
#     plt.plot(punkty[0, :], punkty[1, :], '.', markersize=1)
#     plt.title(rekord.id)
#     plt.show()

# === KONFIGURACJA ===
folder_wejscie = r"C:\Users\Karola\OneDrive - MULTIKKA Sp. z o.o\Dokumenty\STUDIA\INS\magisterka\fasta"
folder_wyjscie = r"C:\Users\Karola\OneDrive - MULTIKKA Sp. z o.o\Dokumenty\STUDIA\INS\magisterka\CGRpng"

# Utwórz folder wyjściowy, jeśli nie istnieje
os.makedirs(folder_wyjscie, exist_ok=True)

# === PRZETWARZANIE WSZYSTKICH FASTA ===
for plik in os.listdir(folder_wejscie):
    if plik.lower().endswith(".fasta"):
        sciezka = os.path.join(folder_wejscie, plik)
        for rekord in SeqIO.parse(sciezka, "fasta"):
            sekwencja = str(rekord.seq)
            punkty = generuj_cgr(sekwencja)

            # Rysowanie i zapisywanie
            plt.figure(figsize=(5,5))
            plt.plot(punkty[0, :], punkty[1, :], '.', markersize=1, color='black')
            plt.axis('off')
            plt.gca().set_position([0, 0, 1, 1])

            nazwa_png = f"{os.path.splitext(plik)[0]}.png"
            sciezka_png = os.path.join(folder_wyjscie, nazwa_png)
            plt.savefig(sciezka_png, dpi=300, bbox_inches='tight', pad_inches=0)
            plt.close()

print("✅ Gotowe! Wszystkie obrazy CGR zapisano w folderze")

print(folder_wyjscie)
