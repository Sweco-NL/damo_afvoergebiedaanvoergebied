import pandas as pd
from importlib_resources import files, as_file

file_waterschappen_codes = files("damo_afvoergebiedaanvoergebied.assets").joinpath("waterschappen_code.csv")
waterschappen_order_codes = pd.read_csv(file_waterschappen_codes, sep=";")
