import os


def fix_csv_header(file_path, output_folder):
    """Fix the header of a single CSV file."""
    try:
        with open(file_path, 'r') as infile:
            lines = infile.readlines()

        # Fix the header
        if len(lines) > 0:
            header = lines[0]
            # Remove quotes and ensure tabs are used
            fixed_header = header.replace('"', '').strip() + '\n'
            rest_of_file = lines[1:]
        else:
            print(f"File {file_path} is empty. Skipping.")
            return

        # Write the fixed file to the output folder
        os.makedirs(output_folder, exist_ok=True)
        output_file = os.path.join(output_folder, os.path.basename(file_path))
        with open(output_file, 'w') as outfile:
            outfile.write(fixed_header)
            outfile.writelines(rest_of_file)

        print(f"Fixed file saved: {output_file}")
    except Exception as e:
        print(f"Error processing {file_path}: {e}")


def process_folder(folder_path, output_folder):
    """Process all CSV files in the folder to fix header formatting."""
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.csv'):
            file_path = os.path.join(folder_path, file_name)
            print(f"Processing file: {file_name}")
            fix_csv_header(file_path, output_folder)


if __name__ == "__main__":
    input_folder = "samCsvs"  # Replace with the folder containing CSV files
    output_folder = "csvs"  # Replace with the folder to save fixed files

    process_folder(input_folder, output_folder)