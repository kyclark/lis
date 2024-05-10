from typing import List, Union
from pydantic import BaseModel
from fastapi import FastAPI
from starlette.middleware.cors import CORSMiddleware
from redis.commands.search import AsyncSearch


class Block(BaseModel):
    chromosome: str
    genus: str
    species: str


#
# Read configuration for global settings
#
config_file = "./config.ini"
assert os.path.isfile(config_file)
config = ConfigParser(interpolation=None)
config.read(config_file)

app = FastAPI()

origins = [
    "http://localhost:*",
    "*",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/blocks", response_model=List[Block])
async def get_blocks(chromosome: str, target: str) -> List[Block]:
    loop = uvloop.new_event_loop()
    conf = config["DEFAULT"]
    connection = await redis.Redis(
        host=conf["redis_host"],
        port=conf["redis_port"],
        db=conf["redis_db"],
        password=conf["redis_password"],
        decode_responses=True
    )
    # ping to force connection, preventing errors downstream
    await connection.ping()
    chromosome_index = AsyncSearch(
        self.redis_connection, index_name="chromosomeIdx"
    )
    return [Block(chromosome="I", genus="Homo", species="Sapiens")]
