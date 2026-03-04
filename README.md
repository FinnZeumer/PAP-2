# PAP 2

## Papulator To Do (for myself lol)

### Wichtige Berechnungen
    [x] Automatisches generieren der Results Datei
        --> Diese soll automatisiert im richtigen Veruschsverzeichnis anhand der Versuchsnummer generiert werden
        --> Ist die File bereits am existieren, so soll diese gelassen werden
        --> Es wird außerdem eine Back-Up File entstehen, die manuell überschreiben werden muss (Button click), damit berechnete Ergebnisse nicht versehentlich gelöscht werden.

    [x] Runden von Fehlern und Messwerten auf signifikante Stellen nach DIN 

    [] Bestimmen der Fehlerformel einer gegebenen Formel
        --> Setzen von fehlerbehafteten Größen und von Exakten Werten

    [] Automatisches berechnen von Messwerten und dazugehörigem Fehler.
        --> Es soll die Funktion gegeben werden. Diese bestimmt dann automatisiert die Fehlerformel. 
        --> Die Messwerte und ihre Ungenauigkeiten sollen ebenfalls gegeben werden. Wichtig ist, das auch das berechnen von mehreren Messwerten möglich sein wird. Heißt wenn man dieselbe Messung wiederholt durchgeführt hat (bspw. unter verschiedenen Bedinugnen), so sollen alle Messungen berechnet werden. Natürlich soll es auch funktionieren, wenn nur eine Messung berechnet wird. 
        --> Es muss unterschieden werden, zwischen 
            a) Fehlerbehafteten Größen
            b) Exakten Werten
            c) Fehlerbehafteten Konstanten
            --> So soll sicher gestellt werden, dass die Fehlerrechnugn gelingt, für Messreihen, die alle dieselbe Konstante brauchen funktionieren. 

### Einlesen von Dateien mit Messwerten
    [x] Speichern eines DataFrames aus einer CSV Datei
    [] TXT Support
    [] Definition eines Bennenungs Schemas.

### Latex Export
