from google import genai
from google.genai import types
from dotenv import load_dotenv
import os


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

# def extract_section_17_information(iteration):
#     iteration = iteration

#     prompt = """
#     You are a mining extraction and refining expert and are evaluating mining process specifications.
#     In section 17 does a specific section exists that evaluates the ion exchange process in the document?
#     If a section for this does exist return integer 1 otherwise return 0
#     In section 17 does a specific section exists that evaluates the flotation process in the document?
#     If a section for this does exist return integer 1 otherwise return 0
#     In section 17 does a specific section exists that evaluates the roasting process in the document?
#     If a section for this does exist return integer 1 otherwise return 0
#     In section 17 does a specific section exists that evaluates the water leaching process in the document?
#     If a section for this does exist return integer 1 otherwise return 0
#     In section 17 does a specific section exists that evaluates the acid leaching process in the document?
#     If a section for this does exist return integer 1 otherwise return 0
#     In section 17 does a specific section exists that evaluates the ion exchange process in the document?
#     If a section for this does exist return integer 1 otherwise return 0
#     In section 17 does a specific section exists that evaluates the precipitation process in the document?
#     If a section for this does exist return integer 1 otherwise return 0
#     In section 17 does a specific section exists that evaluates the calcination process in the document?
#     If a section for this does exist return integer 1 otherwise return 0
#     """

#     response_schema = {
#         "type": "object",
#         "properties": {
#             "exists": {
#                 "type": "object",
#                 "properties": {
#                     "flotation": {"type": "integer"},
#                     "roasting": {"type": "integer"},
#                     "water leaching": {"type": "integer"},
#                     "acid leaching": {"type": "integer"},
#                     "ion_exchange": {"type": "integer"},
#                     "precipitation": {"type": "integer"},
#                     "calcination": {"type": "integer"},
#                 },
#                 "required": ["ion_exchange", "flotation"]
#             }
#         },
#         "required": ["exists"]
#     }

#     client = genai.Client(
#         vertexai=True,
#         project="project-e4f9d719-a476-4261-95d",
#         location="us-central1"
#     )

#     pdf_part = types.Part.from_uri(
#         file_uri=f"gs://mining_literature_test/section_17_chunk_{iteration}.pdf",
#         mime_type="application/pdf"
#     )

#     response = client.models.generate_content(
#         model="gemini-2.5-pro",
#         contents=[pdf_part, prompt],
#         config={
#             "response_mime_type":"application/json",
#             "response_schema":response_schema
#         }
#     )
    
#     return response.text

for i in range(1,4):
    # response = extract_section_17_information(i)
    response = extract_section_17_information(i)
    print(response)

# answers should be 