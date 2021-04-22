import mymodules
import wos

if __name__ == "__main__":
    conn = wos.WOSMongo("../conf/mongo_config.json")
    print("Result:")
    for row in conn.db.journals.find()[0:10]:
        print(row)
