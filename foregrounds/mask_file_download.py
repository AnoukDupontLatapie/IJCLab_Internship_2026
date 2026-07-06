from pathlib import Path
from urllib.request import urlretrieve

url = (
    "https://irsa.ipac.caltech.edu/data/Planck/release_2/"
    "ancillary-data/masks/HFI_Mask_GalPlane-apo0_2048_R2.00.fits"
)


mask_file = Path(__file__).parent / "HFI_Mask_GalPlane-apo0_2048_R2.00.fits"

if not mask_file.exists():
    print("Downloading mask...")
    urlretrieve(url, mask_file)

print(mask_file)
