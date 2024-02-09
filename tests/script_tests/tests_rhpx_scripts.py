import logging
import os
import unittest

import dggstools.__main__ as main
import tempfile
import shutil
from pathlib import Path
from typer.testing import CliRunner



# Some tests to make sure that the scripts are working, and that their user interface
# is not broken.
# Basic functionality is tested in the unit tests.
class RhpxScriptTestsSpec(unittest.TestCase):
    runner = CliRunner(mix_stderr=False) # capture stdout and stderr separately
    cwd = os.getcwd()
    temp_dir = "/tmp" # This should be changed in every test (see the setUp method)
    try:
        data_dir = Path(os.environ['RHPX_TEST_DATA_DIR'])
    except KeyError:
        if os.path.exists("../test_data"):
            data_dir = Path("../test_data")
        elif os.path.exists("test_data"):
            data_dir = Path("test_data")
        elif os.path.exists("tests/test_data"):
            data_dir = Path("tests/test_data")

    # Configuration for the logs (including those of every library used that produces logs)
    logging.basicConfig(format="[%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s", level=logging.INFO)

    test_datasets= ["Aragón_ETRS89_30N.shp", "Aragón_ETRS89_30N.dbf",
                    "Aragón_ETRS89_30N.prj", "Aragón_ETRS89_30N.shx",
                    "landsat_image_small.tif", "NYC_Noise_Part_UTM_18N_WGS84.gpkg"]

    def setUp(self):
        os.chdir(self.cwd)
        self.temp_dir = self._create_tempdir_and_copy_data()
        os.chdir(self.temp_dir)

    def _create_tempdir_and_copy_data(self):
        temp_dir = tempfile.mkdtemp()
        for dataset in self.test_datasets:
            shutil.copyfile(self.data_dir / dataset,
                            Path(temp_dir) / dataset)
        return temp_dir

    def test_v_to_rhealpix(self):
        # Area dataset
        result = self.runner.invoke(main.app, ["vec-to-rhpx-ras",
                                               "Aragón_ETRS89_30N.shp",
                                               "Aragón_ETRS89_30N.tif",
                                               "6"])
        assert result.exit_code == 0
        assert "OK" in result.stdout

        result = self.runner.invoke(main.app, ["vec-to-rhpx-ras",
                                               "Aragón_ETRS89_30N.shp",
                                               "Aragón_ETRS89_30N_nside2.tif",
                                               "8",
                                               "--rdggs", "2/1/0"])
        assert result.exit_code == 0
        assert "OK" in result.stdout

        # Point dataset. It also works, but rasterizing points as this is not very useful
        # Binning the points would be more useful (TODO)
        result = self.runner.invoke(main.app, ["vec-to-rhpx-ras",
                                               "NYC_Noise_Part_UTM_18N_WGS84.gpkg",
                                               "NYC.tif",
                                               "11"])
        print (result.stdout)
        assert result.exit_code == 0
        assert "OK" in result.stdout



    def test_r_to_rhealpix(self):
        result = self.runner.invoke(main.app, ["ras-to-rhpx-ras",
                                               "landsat_image_small.tif",
                                               "landsat_image_small-rHEALPIX.tif",
                                               "11"])
        assert result.exit_code == 0
        assert "OK" in result.stdout

        result = self.runner.invoke(main.app, ["ras-to-rhpx-ras",
                                               "landsat_image_small.tif",
                                               "landsat_image_small-rHEALPIX2.tif",
                                               "12",
                                               "--rdggs", "2/1/0"])
        assert result.exit_code == 0
        assert "OK" in result.stdout


    def test_v_r_area_error(self):
        # First we create the raster
        result = self.runner.invoke(main.app, ["vec-to-rhpx-ras",
                                               "Aragón_ETRS89_30N.shp",
                                               "Aragón_ETRS89_30N_vrareaerr.tif",
                                               "10"])

        result = self.runner.invoke(main.app, ["vec-ras-area-error",
                                                    "Aragón_ETRS89_30N.shp",
                                                    "Aragón_ETRS89_30N_vrareaerr.tif",
                                                    "--input-crs", "EPSG:25830"])
        assert result.exit_code == 0
        assert "RMSE: 690635" in result.stdout

    def test_rhealpix_to_gpkg(self):
        # We need first a rhealpix raster file
        result = self.runner.invoke(main.app, ["ras-to-rhpx-ras",
                                               "landsat_image_small.tif",
                                               "landsat_image_small-rHEALPIX.tif",
                                               "11"])
        assert result.exit_code == 0
        assert "OK" in result.stdout

        result = self.runner.invoke(main.app, ["ras-rhpx-to-vec-rhpx",
                                               "landsat_image_small-rHEALPIX.tif",
                                               "landsat_image_small-rHEALPIX.gpkg"])
        assert result.exit_code == 0
        assert "OK" in result.stdout

        # We can also test the print metadata command which also helps to validate
        # the previous result
        result = self.runner.invoke(main.app, ["print-vec-rhpx-metadata",
                                               "landsat_image_small-rHEALPIX.gpkg"])
        assert result.exit_code == 0
        assert "'res_idx': 11" in result.stdout
        assert "OK" in result.stdout

    def test_gpkg_to_rhealpix(self):
        # We need first a geopackage
        _ = self.runner.invoke(main.app, ["ras-to-rhpx-ras",
                                               "landsat_image_small.tif",
                                               "landsat_image_small-rHEALPIX-921.tif",
                                               "11"])
        _ = self.runner.invoke(main.app, ["ras-rhpx-to-vec-rhpx",
                                               "landsat_image_small-rHEALPIX-921.tif",
                                               "landsat_image_small-rHEALPIX-921.gpkg"])

        result = self.runner.invoke(main.app, ["vec-rhpx-to-ras-rhpx",
                                               "landsat_image_small-rHEALPIX-921.gpkg",
                                               "landsat_image_small-rHEALPIX-921_gpkg.tif"])
        assert result.exit_code == 0
        assert "OK" in result.stdout
