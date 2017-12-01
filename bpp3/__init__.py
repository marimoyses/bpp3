import pkg_resources
import bpp3.bpp3obj
import bpp3.ddPCR
bpp3.helper.subprocess_cmd("module load blatSrc35")
try:
    __version__ = pkg_resources.get_distribution(__name__).version
except:
    __version__ = 'unknown'
