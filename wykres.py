import pandas as pd
import matplotlib.pyplot as plt


FILE_NAME = "taxonomy_labels_FIXED.csv" 

LABEL_MAP = {
    0: "Kotowate",    
    1: "Naczelne",    
    2: "Gryzonie",    
    3: "Walenie",     
    4: "Słoniowate",  
    5: "Psowate",     
    6: "Torbacze",    
    7: "Mroczkowate", 
}

try:
    df = pd.read_csv(FILE_NAME)

    df = df[df['Label_ID'] != -1]

    # Zliczanie 
    class_counts = df['Label_ID'].value_counts().sort_index()

    class_names = [LABEL_MAP.get(i, f"ID {i}") for i in class_counts.index]

    plt.figure(figsize=(10, 6))

    bars = plt.bar(class_names, class_counts.values, color='skyblue')

    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 0.5, yval, ha='center', va='bottom', fontsize=10)

    plt.title('Rozkład Liczby Sekwencji w Poszczególnych Klasach Taksonomicznych', fontsize=14)
    plt.xlabel('Grupa Taksonomiczna', fontsize=12)
    plt.ylabel('Liczba Sekwencji (Próbek)', fontsize=12)

    plt.xticks(rotation=45, ha='right')
    
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

    print("--- Rozkład Liczby Sekwencji ---")
    print(class_counts)
    print("--------------------------------")

except FileNotFoundError:
    print(f"Błąd: Plik '{FILE_NAME}' nie został znaleziony. Upewnij się, że nazwa pliku jest poprawna.")
except KeyError:
    print(f"Błąd: Kolumna 'Label_ID' nie została znaleziona w pliku '{FILE_NAME}'.")