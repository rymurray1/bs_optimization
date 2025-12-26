import pdfplumber
import pandas as pd
import os



# Open the PDF
pdf_path = "/mnt/user-data/uploads/technical_test_doc.pdf"
output_dir = "/mnt/user-data/outputs"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Track all tables found
table_count = 0
tables_info = []

print("Extracting tables from PDF...")
print(f"Processing: {pdf_path}\n")

with pdfplumber.open(pdf_path) as pdf:
    print(f"Total pages: {len(pdf.pages)}\n")
    
    for page_num, page in enumerate(pdf.pages, 1):
        # Extract tables from this page
        tables = page.extract_tables()
        
        if tables:
            print(f"Page {page_num}: Found {len(tables)} table(s)")
            
            for table_num, table in enumerate(tables, 1):
                if table and len(table) > 0:
                    table_count += 1
                    
                    # Create a DataFrame
                    # Use first row as headers if it looks like a header
                    if len(table) > 1:
                        df = pd.DataFrame(table[1:], columns=table[0])
                    else:
                        df = pd.DataFrame(table)
                    
                    # Clean up the DataFrame - remove completely empty rows/columns
                    df = df.dropna(how='all', axis=0)  # Remove empty rows
                    df = df.dropna(how='all', axis=1)  # Remove empty columns
                    
                    # Generate filename
                    filename = f"table_{table_count}_page_{page_num}.csv"
                    filepath = os.path.join(output_dir, filename)
                    
                    # Save to CSV
                    df.to_csv(filepath, index=False)
                    
                    # Store info
                    tables_info.append({
                        'table_number': table_count,
                        'page': page_num,
                        'filename': filename,
                        'rows': len(df),
                        'columns': len(df.columns)
                    })
                    
                    print(f"  - Table {table_num}: {len(df)} rows × {len(df.columns)} columns → {filename}")

print(f"\n{'='*60}")
print(f"Extraction complete!")
print(f"Total tables extracted: {table_count}")
print(f"Output directory: {output_dir}")
print(f"{'='*60}\n")

# Create a summary file
if tables_info:
    summary_df = pd.DataFrame(tables_info)
    summary_path = os.path.join(output_dir, "extraction_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"Summary saved to: extraction_summary.csv")
    print("\nTables extracted:")
    for info in tables_info:
        print(f"  - {info['filename']}: {info['rows']} rows × {info['columns']} columns (Page {info['page']})")
else:
    print("No tables were found in the PDF.")

