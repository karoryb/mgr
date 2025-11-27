import os
import random
from Bio import SeqIO, Seq
from Bio.Data import CodonTable
import pandas as pd

# --- ZMODYFIKOWANE USTAWIANIA FOLDERÓW ---
# Folder, w którym znajdują się pliki FASTA
folder_wejsciowy_fasta = r"C:\Users\Karola\OneDrive - MULTIKKA Sp. z o.o\Dokumenty\STUDIA\INS\magisterka\fasta"
# Folder, w którym znajdują się pliki GFF
folder_wejsciowy_gff = r"C:\Users\Karola\OneDrive - MULTIKKA Sp. z o.o\Dokumenty\STUDIA\INS\magisterka\GFF"
# Folder do zapisania wygenerowanych, zmutowanych sekwencji (pozostaje bez zmian)
folder_wyjsciowy_fasta = r"C:\Users\Karola\OneDrive - MULTIKKA Sp. z o.o\Dokumenty\STUDIA\INS\magisterka\ZMUTOWANE_SEKWENCJE"
# Upewnij się, że folder wyjściowy istnieje
os.makedirs(folder_wyjsciowy_fasta, exist_ok=True)

# Liczba wariantów do wygenerowania na jedną oryginalną sekwencję
NUM_VARIANTS = 10 

# Tabela kodu genetycznego (dla genomu mitochondrialnego kręgowców)
MITO_TABLE = CodonTable.unambiguous_dna_by_id[2]
NUCLEOTIDES = ['A', 'T', 'C', 'G']

# --- WSPÓŁCZYNNIKI MUTACJI ---
# Prawdopodobieństwa mutacji na pozycję (lub kodon)
PROB_CR = 0.005      # Region Kontrolny/D-loop (wysoka zmienność)
PROB_SYN = 0.001     # Mutacja Synonimiczna na kodon (średnia)
PROB_NONSYN = 0.0001 # Mutacja Niesynonimiczna na kodon (niska)
PROB_tRNA_rRNA = 0.00005 # tRNA/rRNA (bardzo niska)
PROB_SNV_OTHER = 0.0005 # Inne regiony nieadnotowane (umiarkowana)

# Bias mutacyjny: Preferowanie przejść (A<->G, C<->T) nad transwersjami
TRANSITION_BIAS = 0.8 # 80% szans na przejście
# Mapowanie przejść
TRANSITIONS = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
# Mapowanie transwersji (dwie możliwości)
TRANSVERSIONS = {'A': ['C', 'T'], 'G': ['C', 'T'], 'C': ['A', 'G'], 'T': ['A', 'G']}

def choose_new_nucleotide(original_nt):
    """Wybiera nowy nukleotyd z uwzględnieniem biasu przejścia/transwersji."""
    if random.random() < TRANSITION_BIAS:
        # Przejście (transition)
        return TRANSITIONS.get(original_nt, random.choice([nt for nt in NUCLEOTIDES if nt != original_nt]))
    else:
        # Transwersja (transversion)
        return random.choice(TRANSVERSIONS.get(original_nt, [nt for nt in NUCLEOTIDES if nt != original_nt]))

def get_synonymous_change(codon):
    """Znajduje losową synonimiczną zmianę dla kodonu, preferując 1 SNV."""
    current_aa = MITO_TABLE.forward_table.get(codon)
    if not current_aa:
        return None  # Kodon STOP lub nieznany

    synonyms = []
    for c, aa in MITO_TABLE.forward_table.items():
        if aa == current_aa and c != codon:
            synonyms.append(c)
    
    # Preferujemy synonimy różniące się tylko 1 nukleotydem (mniej inwazyjne)
    one_snv_synonyms = [s for s in synonyms if sum(s[i] != codon[i] for i in range(3)) == 1]
    
    if one_snv_synonyms:
        return random.choice(one_snv_synonyms)
    elif synonyms:
        return random.choice(synonyms)
    else:
        return None

def parse_gff(gff_file_path):
    """Odczytuje tylko najważniejsze cechy (typy) z pliku GFF za pomocą Pandas."""
    try:
        # Wczytanie pliku GFF, oddzielonego tabulatorami
        df = pd.read_csv(gff_file_path, sep='\t', comment='#', header=None, 
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    except Exception as e:
        print(f"Błąd podczas wczytywania GFF {gff_file_path}: {e}")
        return []

    features = []
    # Interesują nas tylko geny, rRNA, tRNA i D-loop/Control region
    relevant_features = ['gene', 'rRNA', 'tRNA', 'D-loop', 'Control region']
    
    for index, row in df.iterrows():
        feature_type = row['feature']
        if feature_type in relevant_features:
            features.append({
                'type': feature_type,
                # GFF jest 1-bazowe i inkluzywne. Python jest 0-bazowy i ekskluzywny na końcu.
                'start': row['start'] - 1, 
                'end': row['end'], 
                'strand': row['strand']
            })
    return features



def generate_variants(original_seq, features, num_variants):
    """Generuje zmutowane warianty na podstawie sekwencji i adnotacji."""
    variants = []
    seq_len = len(original_seq)
    
    for v_index in range(num_variants):
        mutated_seq_list = list(str(original_seq)) # Kopia do modyfikacji
        
        # 1. Pętla przez całą sekwencję (dla mutacji w regionach nieadnotowanych)
        for pos in range(seq_len):
            # Sprawdzenie, czy pozycja należy do jakiegoś z adnotowanych regionów
            is_annotated = False
            
            # W tym przejściu mutujemy tylko, jeśli pozycja NIE JEST adnotowana 
            # (traktujemy jako inny region niekodujący)
            for feature in features:
                if feature['start'] <= pos < feature['end']:
                    is_annotated = True
                    break
            
            if not is_annotated:
                 if random.random() < PROB_SNV_OTHER:
                    mutated_seq_list[pos] = choose_new_nucleotide(original_seq[pos])

        # 2. Pętla przez adnotowane cechy (zastosowanie precyzyjnych prob.)
        for feature in features:
            f_start = feature['start'] 
            f_end = feature['end'] 
            f_type = feature['type']
            f_strand = feature['strand']

            # Mutacje w regionach niekodujących (D-loop) lub strukturalnych (rRNA/tRNA)
            if f_type in ['D-loop', 'Control region', 'rRNA', 'tRNA']:
                prob = PROB_CR if f_type in ['D-loop', 'Control region'] else PROB_tRNA_rRNA
                for pos in range(f_start, f_end):
                    if random.random() < prob:
                        mutated_seq_list[pos] = choose_new_nucleotide(original_seq[pos])
            
            # Mutacje w regionach kodujących białka ('gene')
            elif f_type == 'gene':
                for codon_start in range(f_start, f_end, 3):
                    codon_end = codon_start + 3
                    if codon_end > f_end:
                        continue # Niepełny kodon na końcu genu
                        
                    # Odczyt kodonu
                    original_codon = "".join(original_seq[codon_start:codon_end])
                    
                    # 1. Mutacja Synonimiczna
                    if random.random() < PROB_SYN:
                        # W przypadku nici minusowej (-), użyj sekwencji komplementarnej
                        if f_strand == '-':
                            temp_codon = str(Seq.Seq(original_codon).reverse_complement())
                        else:
                            temp_codon = original_codon
                            
                        new_codon = get_synonymous_change(temp_codon)
                        
                        if new_codon:
                            # Jeśli gen był na nici minusowej, musimy skomplementować nowo zmutowany kodon
                            if f_strand == '-':
                                final_codon_change = str(Seq.Seq(new_codon).reverse_complement())
                            else:
                                final_codon_change = new_codon

                            # Wprowadzenie synonimicznej zmiany do listy sekwencji
                            for i in range(3):
                                mutated_seq_list[codon_start + i] = final_codon_change[i]
                            continue # Jeśli zaszła synonimiczna, nie próbujemy niesynonimicznej

                    # 2. Mutacja Niesynonimiczna (rzadsza)
                    if random.random() < PROB_NONSYN:
                        # Losowa zmiana nukleotydu w obrębie kodonu
                        nt_to_change_index = random.randint(0, 2) 
                        change_pos = codon_start + nt_to_change_index
                        
                        original_nt = mutated_seq_list[change_pos]
                        mutated_seq_list[change_pos] = choose_new_nucleotide(original_nt)

        variants.append(Seq.Seq("".join(mutated_seq_list)))

    return variants


def main():
    """Główna funkcja przetwarzająca pliki z dwóch różnych folderów."""
    print(f"--- START PRZETWARZANIA ({NUM_VARIANTS} wariantów na gatunek) ---")
    
    # 1. Iterujemy po plikach FASTA (w pełnej nazwie, np. NC_001794.1.fasta)
    fasta_files = [f for f in os.listdir(folder_wejsciowy_fasta) if f.endswith('.fasta') or f.endswith('.fa')]
    
    for fasta_filename in fasta_files:
        # base_name to pełny identyfikator z wersją, np. NC_001794.1
        base_name_full = os.path.splitext(fasta_filename)[0]
        
        # Zidentyfikator bez wersji, np. NC_001794. Użyjemy go do szukania GFF.
        base_name_no_version = base_name_full.split('.')[0] 
        
        fasta_path = os.path.join(folder_wejsciowy_fasta, fasta_filename)
        
        # 2. Szukanie pliku GFF. Oczekujemy nazwy BEZ numeru wersji, np. NC_001794.gff3
        
        # Próba 1: Nazwa bez wersji z .gff3
        gff_filename_gff3 = base_name_no_version + '.gff3'
        gff_path = os.path.join(folder_wejsciowy_gff, gff_filename_gff3)
        
        if not os.path.exists(gff_path):
            # Próba 2: Nazwa bez wersji z .gff
            gff_filename_gff = base_name_no_version + '.gff'
            gff_path = os.path.join(folder_wejsciowy_gff, gff_filename_gff)
            
            if not os.path.exists(gff_path):
                print(f"⚠️ Nie znaleziono pliku GFF dla {base_name_full} (szukano {gff_filename_gff3} lub {gff_filename_gff}). Pomijam.")
                continue
            
        print(f"\n🚀 Przetwarzanie: {base_name_full} (znaleziono GFF: {os.path.basename(gff_path)})")
        
        # 3. Wczytanie i przetwarzanie (reszta kodu bez zmian)
        try:
            record = SeqIO.read(fasta_path, "fasta")
            original_seq = record.seq.upper()
        except Exception as e:
            print(f"❌ Błąd wczytywania FASTA: {e}")
            continue

        features = parse_gff(gff_path)
        if not features:
            print(f"❌ Nie udało się wczytać cech z GFF. Pomijam.")
            continue 

        # Generowanie zmutowanych wariantów
        mutated_variants = generate_variants(original_seq, features, NUM_VARIANTS)
        
        # Zapis wyników (KAŻDY WARIANT OSOBNO)
        save_count = 0
        for i, seq in enumerate(mutated_variants):
            new_id = f"{record.id}_VAR{i+1}"
            new_description = f"Mutowany wariant {i+1} z {record.id}"
            output_record = SeqIO.SeqRecord(seq, id=new_id, description=new_description)
            
            # Nazwa pliku wyjściowego nadal używa pełnego identyfikatora FASTA
            output_fasta_path = os.path.join(folder_wyjsciowy_fasta, f"{base_name_full}_VAR{i+1}.fasta")
            
            try:
                with open(output_fasta_path, "w") as output_handle:
                    SeqIO.write(output_record, output_handle, "fasta")
                save_count += 1
            except Exception as e:
                print(f"❌ Błąd zapisu pliku {output_fasta_path}: {e}")
                
        print(f"✅ Pomyślnie zapisano {save_count} osobnych wariantów dla {base_name_full}.")



if __name__ == "__main__":
    main()