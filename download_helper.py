from urllib.request import urlopen
from io import BytesIO
from zipfile import ZipFile

def download_and_extract_zip(url, extract_to='.'):
    hppt_response = urlopen(url)
    zip_file = ZipFile(BytesIO(hppt_response.read()))
    zip_file.extractall(path=extract_to)

product_group_id = '60835895'
url = 'https://mast.stsci.edu/api/v0.1/Download/bundle.zip?previews=false&obsid=' + product_group_id
destination = '/home/cpp/Área de Trabalho/luz_p/probabilidade/stars'
download_and_extract_zip(url, destination)
