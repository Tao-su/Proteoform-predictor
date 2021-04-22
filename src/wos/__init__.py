import json
from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError

class WOSMongo:
    def __init__(self, filename):
        with open(filename, "rt") as f:
            config = json.loads(f.read())

        uri = WOSMongo.construct_uri(config)
        self.conn = MongoClient(uri)
        self.db = self.conn.web_of_science_aux

    @staticmethod
    def construct_uri(config):
        result = "mongodb://"
        user = config.get("username", None)
        if user is not None:
            password = config.get("password", None)
            result += user
            if password is not None:
                result += ":{}".format(password)
            result += "@"

        host = config["host"]
        port = config.get("port", 27017)
        result += "{}:{}".format(host, port)
            
        result += "/"
        mechanism = config.get("mechanism", None)
        if mechanism is not None:
            if mechanism not in ["MONGODB-CR", "SCRAM-SHA-1"]:
                raise ValueError("mechanism not valid {}".format(mechanism))
            result += "/?authMechanism={}".format(mechanism)
        return result