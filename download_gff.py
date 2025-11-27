import os
import requests
import time
from tqdm import tqdm

# --- USTAWIENIA SKRYPTU ---
EMAIL = "288847@student.pwr.edu.pl"
API_KEY = None # Opcjonalny klucz API NCBI
MAX_TRIES = 5
WAIT_BETWEEN = 0.34 # Kontrola szybkości, aby nie przekroczyć limitów NCBI
INPUT_FILE = "accessions.txt"
OUTPUT_DIR = "GFF" # Zmieniony folder wyjściowy na GFF

# Upewnij się, że folder wyjściowy istnieje
os.makedirs(OUTPUT_DIR, exist_ok=True)


def fetch_gff(acc):
    """Pobiera plik adnotacji GFF3 dla danego accession z NCBI E-utilities."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": acc,
        "rettype": "gff3",  # KLUCZOWA ZMIANA: format GFF3
        "retmode": "text",
        "email": EMAIL,
    }
    if API_KEY:
        params["api_key"] = API_KEY

    headers = {"User-Agent": f"gff-downloader/1.0 ({EMAIL})"}

    for attempt in range(1, MAX_TRIES + 1):
        try:
            # Użycie parametru `timeout` dla bezpieczeństwa
            r = requests.get(base_url, params=params, headers=headers, timeout=30) 
            
            # Sprawdzamy, czy odpowiedź nie jest pusta i ma kod sukcesu
            if r.status_code == 200 and len(r.text) > 50:
                 # GFF/GFF3 zaczynają się od wiersza z metadanymi lub komentarzem ('#')
                 if r.text.startswith("##gff-version 3") or r.text.startswith("#"): 
                    return r.text
                 else:
                    # Czasem NCBI zwraca HTML z błędem, nawet z kodem 200
                    print(f"[{acc}] Serwer zwrócił nieoczekiwaną treść. Sprawdzam, czy accession jest poprawne.")
                    r.raise_for_status() # Wymuszenie błędu
                    
            elif r.status_code in (429, 500, 502, 503, 504):
                sleep_time = 2 ** (attempt - 1) + random.uniform(0, 1)
                print(f"[{acc}] Ograniczenie lub błąd serwera ({r.status_code}), ponawiam za {sleep_time:.2f}s...")
                time.sleep(sleep_time)
            else:
                # Wyrzucenie wyjątku dla innych błędów HTTP
                r.raise_for_status()
                
        except requests.exceptions.RequestException as e:
            sleep_time = 2 ** (attempt - 1) + random.uniform(0, 1)
            print(f"[{acc}] Błąd połączenia lub HTTP: {e}, ponawiam próbę za {sleep_time:.2f}s...")
            time.sleep(sleep_time)

    print(f"[{acc}] ❌ Nie udało się pobrać GFF po {MAX_TRIES} próbach.")
    return None


def main():
    """Główna funkcja wczytująca accessions i pobierająca pliki GFF."""
    try:
        with open(INPUT_FILE, 'r') as f:
            # Usuwanie białych znaków i pustych linii
            accessions = [line.strip().split('.')[0] for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Błąd: Nie znaleziono pliku wejściowego '{INPUT_FILE}'. Upewnij się, że istnieje.")
        return

    print(f"Znaleziono {len(accessions)} identyfikatorów do pobrania.\n")

    # Pobieranie GFF dla każdego accession
    for acc in tqdm(accessions, desc="Pobieranie GFF"):
        # Używamy tylko części przed kropką (np. NC_005089 zamiast NC_005089.1) 
        # do tworzenia nazwy pliku, co jest standardową praktyką.
        
        gff_data = fetch_gff(acc)
        
        if gff_data:
            # Zapis każdego GFF do osobnego pliku
            # Używamy rozszerzenia .gff3, aby pasowało do typu pobranego pliku
            output_path = os.path.join(OUTPUT_DIR, f"{acc}.gff3")
            
            try:
                with open(output_path, "w", encoding="utf-8") as f:
                    f.write(gff_data)
                
            except Exception as e:
                print(f"[{acc}] ❌ Błąd zapisu pliku GFF: {e}")
                
        time.sleep(WAIT_BETWEEN)  # Kontrola szybkości (limity NCBI)

    print(f"\n✅ Zapisano pliki GFF/GFF3 w folderze: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()