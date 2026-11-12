import copy
import unittest

from analysis.palette_evaluation import evaluate, ranks


class PaletteEvaluationTest(unittest.TestCase):
    def setUp(self):
        self.groups = {488: ["a", "b", "c"], 560: ["d", "e", "f"], 640: ["g", "h", "i"]}
        self.values = {}
        for dyes in self.groups.values():
            for i, dye in enumerate(dyes):
                self.values["origami", dye, 1] = self.values["origami", dye, 2] = 9 - i * 3
                self.values["npc", dye, 1] = i + 1
                self.values["npc", dye, 2] = self.values["npc", dye, 3] = 9 - i * 3

    def test_holdout_never_changes_policy(self):
        original = evaluate(self.groups, self.values)
        changed = copy.copy(self.values)
        for dyes in self.groups.values():
            for i, dye in enumerate(dyes):
                changed["npc", dye, 3] = i + 1
                changed["origami", dye, 3] = 1000
                changed["npc", dye, 4] = -5
        updated = evaluate(self.groups, changed)
        self.assertEqual(original["weight"], 1)
        self.assertEqual(original["weight"], updated["weight"])
        self.assertEqual(original["tuning"], updated["tuning"])
        self.assertEqual([c["pair"] for c in original["channels"]], [c["pair"] for c in updated["channels"]])
        self.assertNotEqual(original["macroRegretNm"], updated["macroRegretNm"])

    def test_numerical_ties_and_missing_values(self):
        self.assertEqual(ranks({"a":1,"b":1,"c":3}), {"a":1.5,"b":1.5,"c":3})
        tied = {k: 5.0 for k in self.values}
        result = evaluate(self.groups, tied)
        self.assertEqual(result["weight"], 0)
        self.assertEqual(result["channels"][0]["pair"], ["a", "b"])
        self.assertFalse(result["gatePassed"])
        tied["npc", "i", 3] = None
        result = evaluate(self.groups, tied)
        self.assertEqual(result["status"], "unavailable")
        self.assertIn("npc/i/n3", result["reasons"][0])

    def test_macro_average_and_two_channel_gate(self):
        # Weak repeat-two differences leave the pooled NPC ranking dominated by repeat one.
        for dyes in self.groups.values():
            for i,dye in enumerate(dyes):
                self.values["npc",dye,1] = 5 + 10*i
                self.values["npc",dye,2] = 6 - i
        result = evaluate(self.groups, self.values)
        self.assertTrue(result["gatePassed"])
        self.assertEqual(result["improvedChannels"], 3)
        self.assertEqual(result["macroRegretNm"], 0)
        self.assertEqual(result["baselineMacroRegretNm"], 3)
        for dyes in list(self.groups.values())[1:]:
            for dye in dyes:
                self.values["npc",dye,3] = 1
        one = evaluate(self.groups, self.values)
        self.assertLess(one["macroRegretNm"], one["baselineMacroRegretNm"])
        self.assertEqual(one["improvedChannels"], 1)
        self.assertFalse(one["gatePassed"])
