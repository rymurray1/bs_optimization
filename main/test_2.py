import pdfplumber
import re
from pathlib import Path
import PyPDF2
import shutil
import time
from multiprocessing import Process, Queue


def find_sections_in_pdf(pdf_path, min_page=40, result_queue=None):
    """
    Find Section 13 and Section 17 in a PDF document.

    Args:
        pdf_path: Path to the PDF file
        min_page: Minimum page number to start searching (to skip TOC)
        result_queue: Optional queue for multiprocessing timeout support

    Returns:
        Dictionary with section information:
        {
            'section_13': {'start': int, 'end': int, 'total_pages': int} or None,
            'section_17': {'start': int, 'end': int, 'total_pages': int} or None
        }
    """
    # Section 13 variables
    start_13 = None
    end_13 = None
    used_13_fallback = False  # Track if we used fallback pattern for section 13

    # Section 17 variables
    start_17 = None
    end_17 = None
    used_17_fallback = False  # Track if we used fallback pattern for section 17

    # Match patterns for both sections
    # Primary patterns - look for full section headers
    # (?:\.0)? makes the ".0" optional, [.,:;]? allows optional punctuation after number
    section_13_primary = re.compile(r"^13(?:\.0)?[.,:;]?\s+MINERAL\s+PROCESSING", re.IGNORECASE | re.MULTILINE)
    section_14_primary = re.compile(r"^14(?:\.0)?[.,:;]?\s+\w+", re.IGNORECASE | re.MULTILINE)

    section_17_primary = re.compile(r"^17(?:\.0)?[.,:;]?\s+RECOVERY", re.IGNORECASE | re.MULTILINE)
    section_18_primary = re.compile(r"^18(?:\.0)?[.,:;]?\s+\w+", re.IGNORECASE | re.MULTILINE)

    # Fallback patterns - look for subsection 13.1 or 17.1 (when header text is not connected)
    section_13_fallback = re.compile(r"^13\.1\b", re.IGNORECASE | re.MULTILINE)
    section_14_fallback = re.compile(r"^14\.1\b", re.IGNORECASE | re.MULTILINE)

    section_17_fallback = re.compile(r"^17\.1\b", re.IGNORECASE | re.MULTILINE)
    section_18_fallback = re.compile(r"^18\.1\b", re.IGNORECASE | re.MULTILINE)

    with pdfplumber.open(pdf_path) as pdf:
        for i, page in enumerate(pdf.pages):
            text = page.extract_text() or ""

            # Skip TOC pages
            if i < min_page:
                continue

            # Search for Section 13 (primary pattern with full header)
            if start_13 is None and section_13_primary.search(text):
                start_13 = i
                used_13_fallback = False

            # Search for Section 13 (fallback pattern - 13.1, go back 2 pages)
            elif start_13 is None and section_13_fallback.search(text):
                start_13 = max(min_page, i - 2)  # Go back 2 pages, but not before min_page
                used_13_fallback = True

            # Search for Section 14 (end of Section 13)
            # If we used fallback for 13, look for 14.1; otherwise look for standard 14
            elif start_13 is not None and end_13 is None:
                if used_13_fallback and section_14_fallback.search(text):
                    end_13 = i  # Include the page where 14.1 is found
                elif not used_13_fallback and section_14_primary.search(text):
                    end_13 = i - 1

            # Search for Section 17 (primary pattern with full header)
            if start_17 is None and section_17_primary.search(text):
                start_17 = i
                used_17_fallback = False

            # Search for Section 17 (fallback pattern - 17.1, go back 2 pages)
            elif start_17 is None and section_17_fallback.search(text):
                start_17 = max(min_page, i - 2)  # Go back 2 pages, but not before min_page
                used_17_fallback = True

            # Search for Section 18 (end of Section 17)
            # If we used fallback for 17, look for 18.1; otherwise look for standard 18
            elif start_17 is not None and end_17 is None:
                if used_17_fallback and section_18_fallback.search(text):
                    end_17 = i  # Include the page where 18.1 is found
                elif not used_17_fallback and section_18_primary.search(text):
                    end_17 = i - 1

            # Stop early if we've found both sections completely
            if (start_13 is not None and end_13 is not None and
                start_17 is not None and end_17 is not None):
                break

        # Handle cases where sections extend to end of document
        total_pages = len(pdf.pages)
        if start_13 is not None and end_13 is None:
            end_13 = total_pages - 1

        if start_17 is not None and end_17 is None:
            end_17 = total_pages - 1

    # Build results dictionary
    results = {
        'section_13': None,
        'section_17': None
    }

    if start_13 is not None:
        results['section_13'] = {
            'start': start_13,
            'end': end_13,
            'total_pages': end_13 - start_13 + 1
        }

    if start_17 is not None:
        results['section_17'] = {
            'start': start_17,
            'end': end_17,
            'total_pages': end_17 - start_17 + 1
        }

    # If using multiprocessing, put result in queue
    if result_queue is not None:
        result_queue.put(results)

    return results


def extract_section_to_pdf(pdf_path, start_page, end_page, output_path):
    """
    Extract a specific section from a PDF and save to a new file.

    Args:
        pdf_path: Path to the source PDF file
        start_page: First page to extract (0-indexed)
        end_page: Last page to extract (0-indexed, inclusive)
        output_path: Path for the output PDF file
    """
    with open(pdf_path, 'rb') as file:
        pdf_reader = PyPDF2.PdfReader(file)
        pdf_writer = PyPDF2.PdfWriter()

        # Add pages to the new PDF
        for page_num in range(start_page, end_page + 1):
            pdf_writer.add_page(pdf_reader.pages[page_num])

        # Save the new PDF
        with open(output_path, 'wb') as output_file:
            pdf_writer.write(output_file)


def process_folder(folder_path, min_page=40, split_pdfs=True, timeout_seconds=60):
    """
    Process all PDF files in a folder, find sections, and split them.

    Args:
        folder_path: Path to folder containing PDFs
        min_page: Minimum page number to start searching (to skip TOC)
        split_pdfs: If True, split PDFs into separate section files
        timeout_seconds: Maximum time in seconds to process each PDF (default: 60)

    Returns:
        Dictionary mapping PDF filenames to their section results
    """
    folder = Path(folder_path)
    all_results = {}

    # Create folders if they don't exist
    rejected_folder = folder / "rejected_pdfs"
    split_evaluated_folder = folder / "split_evaluated_pdfs"
    original_evaluated_folder = folder / "original_evaluated_pdfs"

    rejected_folder.mkdir(exist_ok=True)
    split_evaluated_folder.mkdir(exist_ok=True)
    original_evaluated_folder.mkdir(exist_ok=True)

    # Find all PDF files in the folder (using set to avoid duplicates on case-insensitive filesystems)
    pdf_files_lower = list(folder.glob("*.pdf"))
    pdf_files_upper = list(folder.glob("*.PDF"))

    # Combine and remove duplicates by converting to set of paths
    all_pdf_files = list(set(pdf_files_lower + pdf_files_upper))

    # Exclude files in rejected_pdfs, split_evaluated_pdfs, original_evaluated_pdfs folders and already processed files (with _section_13 or _section_17 suffix)
    pdf_files = [f for f in all_pdf_files
                 if "rejected_pdfs" not in str(f)
                 and "split_evaluated_pdfs" not in str(f)
                 and "original_evaluated_pdfs" not in str(f)
                 and not f.stem.endswith("_section_13")
                 and not f.stem.endswith("_section_17")]

    if not pdf_files:
        print(f"No PDF files found in {folder_path}")
        return all_results

    # Create a list snapshot to avoid issues with files being moved during iteration
    pdf_files_snapshot = list(pdf_files)
    total_pdfs_to_process = len(pdf_files_snapshot)
    print(f"Found {total_pdfs_to_process} PDF file(s) to process\n")

    # Process each file in the snapshot
    for pdf_file in pdf_files_snapshot:
        # Check if this PDF has already been processed
        # 1. Original exists in original_evaluated_pdfs
        # 2. Original exists in rejected_pdfs
        # 3. Split files exist in split_evaluated_pdfs
        original_evaluated_copy = original_evaluated_folder / pdf_file.name
        rejected_copy = rejected_folder / pdf_file.name
        split_section_13 = split_evaluated_folder / f"{pdf_file.stem}_section_13.pdf"
        split_section_17 = split_evaluated_folder / f"{pdf_file.stem}_section_17.pdf"

        if original_evaluated_copy.exists() or rejected_copy.exists():
            print(f"Skipping: {pdf_file.name} (already processed - found in original_evaluated or rejected)")
            continue

        if split_section_13.exists() or split_section_17.exists():
            print(f"Skipping: {pdf_file.name} (already processed - split files exist)")
            continue

        print(f"Processing: {pdf_file.name}")
        start_time = time.time()

        try:
            # Create a queue for the result
            result_queue = Queue()

            # Create a process to run the search with timeout
            process = Process(target=find_sections_in_pdf, args=(pdf_file, min_page, result_queue))
            process.start()
            process.join(timeout=timeout_seconds)

            # Check if process completed
            if process.is_alive():
                # Timeout occurred
                process.terminate()
                process.join()
                elapsed_time = time.time() - start_time
                print(f"Processing time: {elapsed_time:.2f}s (TIMEOUT)")
                print(f"Processing exceeded {timeout_seconds}s timeout")

                # Move to rejected folder
                if split_pdfs:
                    rejected_path = rejected_folder / pdf_file.name
                    shutil.move(str(pdf_file), str(rejected_path))
                    print(f"  Moved to rejected_pdfs/{pdf_file.name} (timeout)")

                all_results[pdf_file.name] = {'error': f'Timeout after {timeout_seconds}s'}
                print()
                continue

            # Get results from queue
            if not result_queue.empty():
                results = result_queue.get()
            else:
                raise Exception("No results returned from process")

            elapsed_time = time.time() - start_time
            print(f"  ⏱️  Processing time: {elapsed_time:.2f}s")

            all_results[pdf_file.name] = results

            # Print results for this file
            if results['section_13']:
                sec13 = results['section_13']
                print(f"  Section 13: pages {sec13['start']+1}-{sec13['end']+1} ({sec13['total_pages']} pages)")
            else:
                print(f"  Section 13: Not found")

            if results['section_17']:
                sec17 = results['section_17']
                print(f"  Section 17: pages {sec17['start']+1}-{sec17['end']+1} ({sec17['total_pages']} pages)")
            else:
                print(f"  Section 17: Not found")

            # Split PDFs if both sections found
            if split_pdfs:
                if results['section_13'] and results['section_17']:
                    # Both sections found - split into separate PDFs in split_evaluated folder
                    base_name = pdf_file.stem  # filename without extension

                    # Extract Section 13 to split_evaluated folder
                    section_13_output = split_evaluated_folder / f"{base_name}_section_13.pdf"
                    extract_section_to_pdf(
                        pdf_file,
                        results['section_13']['start'],
                        results['section_13']['end'],
                        section_13_output
                    )
                    print(f"  Created: split_evaluated_pdfs/{section_13_output.name}")

                    # Extract Section 17 to split_evaluated folder
                    section_17_output = split_evaluated_folder / f"{base_name}_section_17.pdf"
                    extract_section_to_pdf(
                        pdf_file,
                        results['section_17']['start'],
                        results['section_17']['end'],
                        section_17_output
                    )
                    print(f"  Created: split_evaluated_pdfs/{section_17_output.name}")

                    # Move original PDF to original_evaluated folder after successful split
                    original_evaluated_path = original_evaluated_folder / pdf_file.name
                    shutil.move(str(pdf_file), str(original_evaluated_path))
                    print(f"  Moved original to original_evaluated_pdfs/{pdf_file.name}")

                else:
                    # One or both sections not found - move to rejected folder
                    rejected_path = rejected_folder / pdf_file.name
                    shutil.move(str(pdf_file), str(rejected_path))
                    print(f"  Moved to rejected_pdfs/{pdf_file.name} (sections not found)")

            print()

        except Exception as e:
            print(f"  Error processing {pdf_file.name}: {e}\n")
            all_results[pdf_file.name] = {'error': str(e)}

            # Move error files to rejected folder
            if split_pdfs:
                rejected_path = rejected_folder / pdf_file.name
                shutil.move(str(pdf_file), str(rejected_path))
                print(f"  Moved to rejected_pdfs/{pdf_file.name} (error occurred)")

    return all_results


if __name__ == "__main__":
    # Example usage: process current directory
    folder_location = "./copper/"
    results = process_folder(folder_location, min_page=40, split_pdfs=True)

    print("\n=== SUMMARY ===")
    successful = 0
    rejected = 0

    for filename, data in results.items():
        print(f"\n{filename}:")
        if 'error' in data:
            print(f"  Error: {data['error']}")
            rejected += 1
        else:
            if data['section_13'] and data['section_17']:
                print(f"  Section 13: Found")
                print(f"  Section 17: Found")
                print(f"  Status: Successfully split into separate PDFs")
                successful += 1
            else:
                if data['section_13']:
                    print(f"  Section 13: Found")
                else:
                    print(f"  Section 13: Not found")

                if data['section_17']:
                    print(f"  Section 17: Found")
                else:
                    print(f"  Section 17: Not found")

                print(f"  Status: Moved to rejected_pdfs/")
                rejected += 1

    print(f"\n=== FINAL STATS ===")
    print(f"Successfully processed: {successful}")
    print(f"Rejected: {rejected}")
    print(f"Total: {successful + rejected}")