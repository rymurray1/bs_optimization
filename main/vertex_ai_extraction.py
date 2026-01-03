from google import genai
from google.genai import types
from dotenv import load_dotenv
import os
import json
import csv
import time


load_dotenv()
api_key = os.getenv("GEMINI_API_KEY")


def extract_section_13_information(iteration):
    iteration = iteration

    prompt = """
    You are a mining extraction and refining expert.
    Evaluate section 13 and return the copper recovery percentage and the grind size.
    Do not evaluate a section other than 13.
    If either value cannot be found, return null.
    Your response should be limited to the most relevant number only. 
    Return only the required value and nothing more. 
    If a previous response contains an answer that is more accurate than what you discovered, prefer the previous response and return null.
    """

    response_schema = {
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "copper_recovery": {"type": "number"},
                "grind_size": {"type": "integer"},
            },
            "required": ["copper_recovery", "grind_size"],
        },
    }

    client = genai.Client(
        vertexai=True,
        project="project-e4f9d719-a476-4261-95d",
        location="us-central1"
    )

    pdf_part = types.Part.from_uri(
        file_uri=f"gs://mining_literature_test/chunk_{iteration}.pdf",
        mime_type="application/pdf"
    )

    response = client.models.generate_content(
        model="gemini-3-pro",
        contents=[pdf_part, prompt],
        config={
            "response_mime_type":"application/json",
            "response_schema":response_schema
        }
    )
    
    return response.text


def extract_section_17_information(iteration):
    iteration = iteration

    prompt = """
    You are a mining extraction and refining expert and are evaluating mining process specifications.
    Your goal is to find and return information for each stage of importance. 
    You should only evaluate section 17 of the document under Recovery Methods.
    Do not evaluate a section other than 17.
    The important stages are: grinding/crushing, flotation, roasting, water leaching, acid leaching, ion exchange, precipitation, calcination.
    Is grinding/crushing discussed in the document?
    If so return the abrasiveness, operating speed, reduction ratio and grind size else return null.
    Return only the information for grinding/crushing you are certain of.
    Is flotation discussed in the document?
    If so return pulp density, reagents, frother type, airflow rate, feed flowrate, wash water, flowrate and retention time else return null.
    Return only the information for flotation you are certain of. 
    Is roasting discussed in the document?
    If so return solid to acid ratio, temperature, and time else return null.
    Return only the information for water leaching you are certain of.  
    Is water leaching discussed in the document?
    If so return pulp density, temperature and time else return null.
    Is acid leaching discussed in the document?
    If so return pulp density, acid, temperature and time else return null
    Is ion exchange discussed in the document?
    If so return pH and loading else return null.
    Return only the information for precipitation you are certain of.
    Is precipitation discussed in the document?
    If so return acid/feed as a ratio, temperature, time and pressure else return null.
    Return only the information for roasting you are certain of.
    Is calcination discussed in the document?
    If so return temperature and pressure.
    Return null if you aren't able to find a suitable answer.
    If any number is evaluated as a range, return the halfway point between the two numbers. 
    Your response should be limited to the most relevant number or string for each section only. 
    Return only the required value and nothing more. 
    If a previous response was more relevant, return null.
    """

    response_schema = {
        "type": "object",
        "properties": {
            "grinding_crushing": {
                "type": "object",
                "properties": {
                    "abrasiveness": {"type": "number"},
                    "operating_speed": {"type": "number"},
                    "reduction_ratio": {"type": "number"},
                    "grind_size": {"type": "number"},
                },
                "required": ["abrasiveness", "operating_speed", "reduction_ratio", "grind_size"],
            },
            "flotation": {
                "type": "object",
                "properties": {
                    "pulp_density": {"type": "integer"},
                    "reagents": {"type": "string"},
                    "frother": {"type": "string"},
                    "airflow_rate": {"type": "number"},
                    "wash_water": {"type": "number"},
                    "flowrate": {"type": "number"},
                    "retention_time": {"type": "integer"},
                },
                "required": ["pulp_density", "reagents", "frother", "airflow_rate", "wash_water", "flowrate", "retention_time"]
            },
            "roasting": {
                "type": "object",
                "properties": {
                    "acid_to_solid_ratio": {"type": "number"},
                    "temperature": {"type": "number"},
                    "time": {"type": "integer"},
                },
                "required": ["acid_to_solid_ratio", "temperature", "time"]
            },
            "water_leaching": {
                "type": "object",
                "properties": {
                    "pulp_density": {"type": "integer"},
                    "temperature": {"type": "integer"},
                    "time": {"type": "number"},
                },
                "required": ["pulp_density", "temperature", "time"]
            },
            "acid_leaching": {
                "type": "object",
                "properties": {
                    "pulp_density": {"type": "integer"},
                    "acid": {"type": "string"},
                    "temperature": {"type": "integer"},
                    "time": {"type": "number"},
                },
                "required": ["pulp_density", "acid", "temperature", "time"]
            },
            "ion_exchange": {
                "type": "object",
                "properties": {
                    "ph": {"type": "number"},
                    "loading": {"type": "number"},
                },
                "required": ["ph", "loading"]
            },
            "precipitation": {
                "type": "object",
                "properties": {
                    "acid_to_feed_ratio": {"type": "number"},
                    "temperature": {"type": "integer"},
                    "time": {"type": "number"},
                    "pressure": {"type": "string"},
                },
                "required": ["acid_to_feed_ratio", "temperature", "time", "pressure"]
            },
            "calcination": {
                "type": "object",
                "properties": {
                    "temperature": {"type": "number"},
                    "pressure": {"type": "string"},
                },
                "required": ["temperature", "pressure"]
            }
        },
        "required": ["grinding_crushing", "flotation", "roasting", "water_leaching", "acid_leaching", "ion_exchange", "precipitation", "calcination"]
    }

    client = genai.Client(
        vertexai=True,
        project="project-e4f9d719-a476-4261-95d",
        location="us-central1"
    )

    pdf_part = types.Part.from_uri(
        file_uri=f"gs://mining_literature_test/chunk_{iteration}.pdf",
        mime_type="application/pdf"
    )

    response = client.models.generate_content(
        model="gemini-2.5-pro",
        contents=[pdf_part, prompt],
        config={
            "response_mime_type":"application/json",
            "response_schema":response_schema
        }
    )
    
    return response.text


def extract_mineral_processing_data(gcs_uri):
    """
    Extract mineral processing data from Section 13 or Section 17.

    Args:
        gcs_uri: GCS URI to the PDF file (e.g., "gs://mining_literature_test/document.pdf")

    Returns:
        JSON string with extracted data
    """
    # Add unique timestamp to prevent caching
    import uuid
    unique_id = str(uuid.uuid4())

    prompt = f"""
    You are a mineral processing and mining data extraction expert.
    Extract mineral processing operational parameters ONLY from Section 17 (Recovery Methods) of THIS SPECIFIC document.

    IMPORTANT: Do NOT use any information from previous documents or previous responses. Only extract data that is explicitly present in THIS document.
    Each document is unique and must be evaluated independently.

    Request ID: {unique_id}
    Document: {gcs_uri}

    For each parameter below, search for the value and units in the document:
    - throughput (units: t/d)
    - head_grade (units: %)
    - crushing_availability (units: %)
    - mill_availability (units: %)
    - concentrate_availability (units: %)
    - tailings_filtration_availability (units: %)
    - bond_crusher_work_index (units: metric)
    - bond_mill_work_index (units: metric)
    - smc_axb (dimensionless)
    - bond_abrasion_index (units: g)
    - ore_specific_gravity (dimensionless)
    - sag_mill_discharge_density (units: % w/w)
    - sag_mill_discharge_availability (units: % v/v)
    - pebble_crusher_recycle_rate (units: % w/w)
    - ball_mill_discharge_rate (units: % w/w)
    - ball_mill_discharge_density (units: % v/v)
    - gravity_circuit_feed_rate (units: %)
    - grind_size (units: µm)
    - class_circulating_load (units: %)
    - rougher_stage_recovery (units: %)
    - rougher_stage_mass_recovery (units: %)
    - cleaner_stage_recovery (units: %)
    - cleaner_stage_mass_recovery (units: %)
    - regrind_size (units: µm)
    - concentrate_thickener_settling_rate (units: t/m2/h)
    - concentrate_thickener_underflow_density (units: % w/w)
    - tailings_thickener_settling_rate (units: t/m2/h)
    - tailing_thickener_underflow_density (units: % w/w)
    - concentrate_moisture_content (units: % w/w)
    - tailings_filter_cake_moisture_content (units: % w/w)
    - sulphide_tails_filter_cake_moisture_content (units: % w/w)
    - feed_size (variable units)
    - frother_consumption (units: t/y)
    - reagent_consumption (variable units)
    - reagent_used (text/name)
    - frother_used (text/name)
    - cyclone_overflow_solids (variable units)
    - flotation_feed_density (variable units)
    - metal_recovery_copper (units: %)
    - metal_concentrate_copper (units: %)

    Return null for any string parameter not found in the document and -1 for any number not found.
    If a value is found, extract the numeric value only (for numeric fields) or the text/name (for text fields).
    If a range is given, return the midpoint.
    Be precise and only return values you are confident about from the document.
    """

    response_schema = {
        "type": "object",
        "properties": {
            "throughput": {"type": "number"},
            "head_grade": {"type": "number"},
            "crushing_availability": {"type": "number"},
            "mill_availability": {"type": "number"},
            "concentrate_availability": {"type": "number"},
            "tailings_filtration_availability": {"type": "number"},
            "bond_crusher_work_index": {"type": "number"},
            "bond_mill_work_index": {"type": "number"},
            "smc_axb": {"type": "number"},
            "bond_abrasion_index": {"type": "number"},
            "ore_specific_gravity": {"type": "number"},
            "sag_mill_discharge_density": {"type": "number"},
            "sag_mill_discharge_availability": {"type": "number"},
            "pebble_crusher_recycle_rate": {"type": "number"},
            "ball_mill_discharge_rate": {"type": "number"},
            "ball_mill_discharge_density": {"type": "number"},
            "gravity_circuit_feed_rate": {"type": "number"},
            "grind_size": {"type": "number"},
            "class_circulating_load": {"type": "number"},
            "rougher_stage_recovery": {"type": "number"},
            "rougher_stage_mass_recovery": {"type": "number"},
            "cleaner_stage_recovery": {"type": "number"},
            "cleaner_stage_mass_recovery": {"type": "number"},
            "regrind_size": {"type": "number"},
            "concentrate_thickener_settling_rate": {"type": "number"},
            "concentrate_thickener_underflow_density": {"type": "number"},
            "tailings_thickener_settling_rate": {"type": "number"},
            "tailing_thickener_underflow_density": {"type": "number"},
            "concentrate_moisture_content": {"type": "number"},
            "tailings_filter_cake_moisture_content": {"type": "number"},
            "sulphide_tails_filter_cake_moisture_content": {"type": "number"},
            "feed_size": {"type": "string"},
            "frother_consumption": {"type": "number"},
            "reagent_consumption": {"type": "string"},
            "reagent_used": {"type": "string"},
            "frother_used": {"type": "string"},
            "cyclone_overflow_solids": {"type": "string"},
            "flotation_feed_density": {"type": "string"},
            "metal_recovery_copper": {"type": "number"},
            "metal_concentrate_copper": {"type": "number"}
        }
    }

    client = genai.Client(
        vertexai=True,
        project="project-e4f9d719-a476-4261-95d",
        location="us-central1"
    )

    pdf_part = types.Part.from_uri(
        file_uri=gcs_uri,
        mime_type="application/pdf"
    )

    response = client.models.generate_content(
        model="gemini-2.5-pro",
        contents=[pdf_part, prompt],
        config={
            "response_mime_type": "application/json",
            "response_schema": response_schema,
            "temperature": 0.0,
            "top_p": 1.0,
        }
    )

    # Get the result before closing the client
    result = response.text

    # Close the client to ensure no session persistence between extractions
    del client

    return result


def check_file_already_processed(output_file, filename):
    """
    Check if a file has already been processed by looking in the CSV.

    Args:
        output_file: Path to the CSV file
        filename: Name of the file to check

    Returns:
        True if file already exists in CSV, False otherwise
    """
    if not os.path.exists(output_file):
        return False

    try:
        with open(output_file, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if row.get('file_name') == filename:
                    return True
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return False

    return False


def export_to_csv(json_data, output_file="vertex_extraction.csv", id_number=101, filename=""):
    """
    Export JSON data to CSV file.

    Args:
        json_data: JSON string or dictionary containing the extracted data
        output_file: Name of the output CSV file (default: vertex_extraction.csv)
        id_number: ID number to identify the data source (default: 101)
        filename: Name of the source file
    """
    # Parse JSON if it's a string
    if isinstance(json_data, str):
        data = json.loads(json_data)
    else:
        data = json_data

    # Add ID and file_name to the data
    data_with_id = {"id": id_number, "file_name": filename}
    data_with_id.update(data)

    # Check if file exists and has data
    file_exists = os.path.exists(output_file)

    # Write to CSV
    with open(output_file, 'a' if file_exists else 'w', newline='', encoding='utf-8') as csvfile:
        # Get all field names from the data (id and file_name will be first)
        fieldnames = list(data_with_id.keys())

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write header only if file doesn't exist
        if not file_exists:
            writer.writeheader()

        # Write data row
        writer.writerow(data_with_id)

    print(f"\nData exported to {output_file}")


# Example usage
if __name__ == "__main__":
    # Use the manually uploaded file

    for i in range(4, 7):
        gcs_uri = f"gs://mining_literature_test/test_document_{i}.pdf"
        filename = f"test_document_{i}.pdf"
        output_file = "vertex_extraction.csv"

        # Check if file has already been processed
        if check_file_already_processed(output_file, filename):
            print(f"\n'{filename}' has already been processed. Skipping extraction.")
        else:
            print(f"Extracting mineral processing data from {gcs_uri}...")
            result = extract_mineral_processing_data(gcs_uri)
            print("\nExtracted Data:")
            print(result)

            # Export to CSV
            export_to_csv(result, output_file, id_number=101, filename=filename)
    # gcs_uri = "gs://mining_literature_test/test_document_4.pdf"
    # filename = "test_document_4.pdf"
    # output_file = "vertex_extraction.csv"

    # # Check if file has already been processed
    # if check_file_already_processed(output_file, filename):
    #     print(f"\n'{filename}' has already been processed. Skipping extraction.")
    # else:
    #     print(f"Extracting mineral processing data from {gcs_uri}...")
    #     result = extract_mineral_processing_data(gcs_uri)
    #     print("\nExtracted Data:")
    #     print(result)

    #     # Export to CSV
    #     export_to_csv(result, output_file, id_number=101, filename=filename)