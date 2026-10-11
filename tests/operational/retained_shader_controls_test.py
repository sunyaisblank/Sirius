"""Keep optional SPIR-V capability names intact when adding float controls."""

import importlib.util
from pathlib import Path
import re
import sys
import unittest

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "retained_shader_builder", ROOT / "scripts" / "build-retained-kernels.py")
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)


class RetainedShaderControlsTest(unittest.TestCase):
    def module(self, capabilities):
        return "\n".join("               OpCapability " + name for name in capabilities) + """
               OpEntryPoint GLCompute %main "main"
               OpExecutionMode %main LocalSize 64 1 1
               OpExecutionMode %main DenormPreserve 32
               OpDecorate %sum NoContraction
       %sum = OpFAdd %float %a %b
"""

    def test_shader_clock_capability_survives_either_declaration_order(self):
        # ShaderClockKHR and Shader are distinct SPIR-V capability tokens.
        for capabilities in (("ShaderClockKHR", "Shader", "DenormPreserve"),
                             ("Shader", "DenormPreserve", "ShaderClockKHR")):
            with self.subTest(capabilities=capabilities):
                result = builder.native_controls(self.module(capabilities))
                actual = re.findall(r"OpCapability (\S+)", result)
                self.assertCountEqual(actual, (*capabilities, "RoundingModeRTE"))
                self.assertIn("OpExecutionMode %main RoundingModeRTE 32", result)
                self.assertIn("OpDecorate %sum NoContraction", result)

    def test_original_native_module_retains_its_required_controls(self):
        result = builder.native_controls(self.module(("Shader", "DenormPreserve")))
        self.assertCountEqual(re.findall(r"OpCapability (\S+)", result),
                              ("Shader", "DenormPreserve", "RoundingModeRTE"))
        self.assertIn("OpExecutionMode %main DenormPreserve 32", result)
        self.assertIn("OpExecutionMode %main RoundingModeRTE 32", result)
        self.assertIn("%sum = OpFAdd %float %a %b", result)


if __name__ == "__main__":
    unittest.main()
