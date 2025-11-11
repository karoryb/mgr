import os
import requests
import time
from tqdm import tqdm


EMAIL = "288847@student.pwr.edu.pl"
API_KEY = None
MAX_TRIES = 5
WAIT_BETWEEN = 0.34
INPUT_FILE = "accessions.txt"
OUTPUT_MODE = "separate" 
COMBINED_FILE = "all_sequences.fasta"
OUTPUT_DIR = "fasta" 


os.makedirs(OUTPUT_DIR, exist_ok=True)


def fetch_fasta(acc):
    """Pobiera sekwencję FASTA dla danego accession z NCBI E-utilities"""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": acc,
        "rettype": "fasta",
        "retmode": "text",
        "email": EMAIL,
    }
    if API_KEY:
        params["api_key"] = API_KEY

    headers = {"User-Agent": f"fasta-downloader/1.0 ({EMAIL})"}

    for attempt in range(1, MAX_TRIES + 1):
        try:
            r = requests.get(base_url, params=params, headers=headers, timeout=20)
            if r.status_code == 200 and r.text.startswith(">"):
                return r.text
            elif r.status_code in (429, 500, 502, 503, 504):
                sleep_time = 2 ** (attempt - 1)
                print(f"[{acc}] Ograniczenie lub błąd serwera ({r.status_code}), ponawiam za {sleep_time}s...")
                time.sleep(sleep_time)
            else:
                r.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"[{acc}] Błąd: {e}, ponawiam próbę...")
            time.sleep(2 ** (attempt - 1))
    print(f"[{acc}] ❌ Nie udało się pobrać po {MAX_TRIES} próbach.")
    return None


def main():
    with open(INPUT_FILE) as f:
        accessions = [line.strip() for line in f if line.strip()]

    print(f"Znaleziono {len(accessions)} identyfikatorów do pobrania.\n")

    if OUTPUT_MODE == "combined":
        combined_path = os.path.join(OUTPUT_DIR, COMBINED_FILE)
        combined = open(combined_path, "w", encoding="utf-8")

    for acc in tqdm(accessions, desc="Pobieranie FASTA"):
        seq = fetch_fasta(acc)
        if seq:
            if OUTPUT_MODE == "combined":
                combined.write(seq + "\n")
            else:
                output_path = os.path.join(OUTPUT_DIR, f"{acc}.fasta")
                with open(output_path, "w", encoding="utf-8") as f:
                    f.write(seq)
        time.sleep(WAIT_BETWEEN)  # kontrola szybkości (limity NCBI)

    if OUTPUT_MODE == "combined":
        combined.close()
        print(f"\nWszystkie sekwencje zapisano do: {combined_path}")
    else:
        print(f"\nSekwencje zapisano w folderze: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
