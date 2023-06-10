import unittest
import pathlib
import os

import gurobipy as gp
from model_analyzer import kappa_explain, angle_explain

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
        kappa_explain(self.model, expltype="ROWS")
        assert lpfile.exists()

    def test_bycols(self):
        mpsfile = cwd.joinpath("testmodel_kappaexplain.mps")
        assert not mpsfile.exists()
        kappa_explain(self.model, expltype="COLS")
        assert mpsfile.exists()

    def test_lpsubprob(self):
        lpfile = cwd.joinpath("testmodel_kappaexplain.lp")
        assert not lpfile.exists()
        kappa_explain(self.model, relobjtype="LP")
        assert lpfile.exists()

    def test_qpsubprob(self):
        lpfile = cwd.joinpath("testmodel_kappaexplain.lp")
        assert not lpfile.exists()
        kappa_explain(self.model, relobjtype="QP")
        assert lpfile.exists()

    def test_default(self):
        lpfile = cwd.joinpath("testmodel_kappaexplain.lp")
        assert not lpfile.exists()
        kappa_explain(self.model, method="DEFAULT")
        assert lpfile.exists()

    def test_lasso(self):
        lpfile = cwd.joinpath("testmodel_kappaexplain.lp")
        assert not lpfile.exists()
        kappa_explain(self.model, method="LASSO")
        assert lpfile.exists()

    def test_rls(self):
        lpfile = cwd.joinpath("testmodel_kappaexplain.lp")
        assert not lpfile.exists()
        kappa_explain(self.model, method="RLS")
        assert lpfile.exists()

    def test_submatrix(self):
        lpfile = cwd.joinpath("testmodel_kappaexplain.lp")
        assert not lpfile.exists()
        kappa_explain(self.model, submatrix=True)
        assert lpfile.exists()
    #
    #   angle_explain is different, as it doesn't write a file
    #   There's a problem is any of the 3 items it returns
    #   is None
    #
    def test_angles(self):
        list1, list2, model = kappa_explain(self.model, method="ANGLES")
        assert list1 != None and list2 != None and model != None
