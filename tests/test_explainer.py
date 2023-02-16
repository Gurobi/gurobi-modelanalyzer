import unittest
import pathlib
import os

import gurobipy as gp
from model_analyzer import kappa_explain, BYROWS, BYCOLS

here = pathlib.Path(__file__).parent
cwd = pathlib.Path(os.getcwd())


class TestExplainer(unittest.TestCase):
    """ Super simple tests: did we write output files? """

    def setUp(self):
        self.clean()
        self.env = gp.Env()
        self.mip_model = gp.read(str(here / "dataset" / "p0033.lp"), env=self.env)
        self.model = self.mip_model.relax()
        self.model.ModelName = "testmodel"

    def tearDown(self):
        self.model.close()
        self.mip_model.close()
        self.env.close()
        self.clean()

    def clean(self):
        output_files = [
            cwd.joinpath("testmodel_kappaexplain.lp"),
            cwd.joinpath("testmodel_kappaexplain.mps"),
        ]
        for file in output_files:
            if file.exists():
                os.remove(file)

    def test_byrows(self):
        lpfile = cwd.joinpath("testmodel_kappaexplain.lp")
        assert not lpfile.exists()
        kappa_explain(self.model, expltype=BYROWS)
        assert lpfile.exists()

    def test_bycols(self):
        mpsfile = cwd.joinpath("testmodel_kappaexplain.mps")
        assert not mpsfile.exists()
        kappa_explain(self.model, expltype=BYCOLS)
        assert mpsfile.exists()
