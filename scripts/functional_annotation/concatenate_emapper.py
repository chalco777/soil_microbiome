import os

# Directorio con los archivos GFF a unir
input_dir = "/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/gff_final_3/gff"
output_file = "/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/gff_final_3/combined_gff_final.gff"

# Crear o sobrescribir el archivo de salida
with open(output_file, 'w') as outfile:
    # Iterar sobre todos los archivos en el directorio
    for filename in sorted(os.listdir(input_dir)):
        if filename.endswith(".gff"):
            barcode_name = filename.split('_')[0]  # Extraer el nombre del barcode
            file_path = os.path.join(input_dir, filename)
            
            # Leer y procesar cada archivo GFF
            with open(file_path, 'r') as infile:
                for line in infile:
                    # Ignorar líneas en blanco
                    if line.strip():
                        # Añadir el nombre del barcode como primera columna
                        if not line.startswith("#"):
                            outfile.write(f"{barcode_name}\t{line}")
                        else:
                            # Mantener comentarios sin añadir la columna del barcode
                            outfile.write(line)

print(f"Archivo GFF combinado creado en: {output_file}")
