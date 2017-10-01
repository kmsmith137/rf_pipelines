from web_viewer import app
from web_viewer import Crawler

if __name__ == "__main__":
    master_directories = Crawler()
    app.run()
