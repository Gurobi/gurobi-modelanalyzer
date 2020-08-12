from unittest import TestCase
import os

import model_analyzer.analyzer as analyzer


class TestAnalyzer(TestCase):
    def setUp(self) -> None:
        cwd = os.getcwd()

        self.file_path = os.path.join(cwd, 'dataset', 'glass4.mps')

    def test_run(self):
        data = analyzer.run(self.file_path)
        print("Done")
