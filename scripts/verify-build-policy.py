#!/usr/bin/env python3
"""Enforce Sirius's non-negotiable qualification and release build topology."""

from __future__ import annotations

import argparse
import sys


def parse_bool(value: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise ValueError(f"expected true or false, got {value!r}")


def configure_errors(values: dict[str, object]) -> list[str]:
    if values["alignment_mode"] not in {"qualification", "release"}:
        return []
    errors = []
    required_true = {
        "build_tests": "BUILD_TESTS",
        "mandatory_tests": "SIRIUS_MANDATORY_TESTS",
        "warnings_as_errors": "SIRIUS_WERROR",
        "require_vulkan_runtime": "SIRIUS_REQUIRE_VULKAN_RUNTIME",
    }
    for key, cmake_name in required_true.items():
        if values[key] is not True:
            errors.append(f"strict policy requires {cmake_name}=ON")
    if values["contract_mode"] != "2":
        errors.append("strict policy requires SIRIUS_CONTRACT_MODE=2")

    configurations = values["configurations"]
    if not isinstance(configurations, tuple) or configurations != ("Release",):
        errors.append("strict policy permits only the Release build configuration")
    return errors


def product_errors(values: dict[str, object]) -> list[str]:
    if values["alignment_mode"] not in {"qualification", "release"}:
        return []
    errors = []
    required_true = {
        "vulkan_backend": "the Vulkan backend",
        "compiled_kernels": "compiled Slang kernels",
        "spirv_validation": "SPIR-V validation",
        "native_viewer": "the native viewer",
    }
    for key, description in required_true.items():
        if values[key] is not True:
            errors.append(f"strict product requires {description}")
    return errors


def self_test() -> None:
    configure = {
        "alignment_mode": "release",
        "build_tests": True,
        "mandatory_tests": True,
        "warnings_as_errors": True,
        "require_vulkan_runtime": True,
        "contract_mode": "2",
        "configurations": ("Release",),
    }
    product = {
        "alignment_mode": "release",
        "vulkan_backend": True,
        "compiled_kernels": True,
        "spirv_validation": True,
        "native_viewer": True,
    }
    if configure_errors(configure) or product_errors(product):
        raise ValueError("complete release policy was rejected")
    qualification_configure = dict(configure, alignment_mode="qualification")
    qualification_product = dict(product, alignment_mode="qualification")
    if configure_errors(qualification_configure) or product_errors(qualification_product):
        raise ValueError("strict evidence-qualification policy was rejected")
    for strict_mode in ("qualification", "release"):
        for key in (
            "build_tests",
            "mandatory_tests",
            "warnings_as_errors",
            "require_vulkan_runtime",
        ):
            weakened = dict(configure, alignment_mode=strict_mode)
            weakened[key] = False
            if not configure_errors(weakened):
                raise ValueError(f"{strict_mode} policy accepted false {key}")
    for contract_mode in ("0", "1"):
        weakened = dict(configure)
        weakened["contract_mode"] = contract_mode
        if not configure_errors(weakened):
            raise ValueError(f"strict policy accepted contract mode {contract_mode}")
    for configurations in ((), ("Debug",), ("Debug", "Release")):
        weakened = dict(configure)
        weakened["configurations"] = configurations
        if not configure_errors(weakened):
            raise ValueError(f"strict policy accepted configurations {configurations}")
    for strict_mode in ("qualification", "release"):
        for key in ("vulkan_backend", "compiled_kernels", "spirv_validation", "native_viewer"):
            weakened = dict(product, alignment_mode=strict_mode)
            weakened[key] = False
            if not product_errors(weakened):
                raise ValueError(f"{strict_mode} product accepted missing {key}")
    development = dict(configure)
    development.update(
        alignment_mode="development",
        build_tests=False,
        mandatory_tests=False,
        warnings_as_errors=False,
        require_vulkan_runtime=False,
        contract_mode="0",
        configurations=("Debug", "Release"),
    )
    if configure_errors(development):
        raise ValueError("development policy was incorrectly promoted to release")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--phase", choices=("configure", "product"))
    parser.add_argument("--alignment-mode", choices=("development", "qualification", "release"))
    parser.add_argument("--build-tests", choices=("true", "false"))
    parser.add_argument("--mandatory-tests", choices=("true", "false"))
    parser.add_argument("--warnings-as-errors", choices=("true", "false"))
    parser.add_argument("--require-vulkan-runtime", choices=("true", "false"))
    parser.add_argument("--contract-mode")
    parser.add_argument("--configurations")
    parser.add_argument("--vulkan-backend", choices=("true", "false"))
    parser.add_argument("--compiled-kernels", choices=("true", "false"))
    parser.add_argument("--spirv-validation", choices=("true", "false"))
    parser.add_argument("--native-viewer", choices=("true", "false"))
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            print("strict qualification/release policy rejected every weakened configuration and product")
            return 0
        if args.phase is None or args.alignment_mode is None:
            raise ValueError("--phase and --alignment-mode are required")
        if args.phase == "configure":
            required = (
                args.build_tests,
                args.mandatory_tests,
                args.warnings_as_errors,
                args.require_vulkan_runtime,
                args.contract_mode,
                args.configurations,
            )
            if any(value is None for value in required):
                raise ValueError("configure phase is missing required arguments")
            values = {
                "alignment_mode": args.alignment_mode,
                "build_tests": parse_bool(args.build_tests),
                "mandatory_tests": parse_bool(args.mandatory_tests),
                "warnings_as_errors": parse_bool(args.warnings_as_errors),
                "require_vulkan_runtime": parse_bool(args.require_vulkan_runtime),
                "contract_mode": args.contract_mode,
                "configurations": tuple(
                    item for item in args.configurations.split(",") if item
                ),
            }
            errors = configure_errors(values)
        else:
            required = (
                args.vulkan_backend,
                args.compiled_kernels,
                args.spirv_validation,
                args.native_viewer,
            )
            if any(value is None for value in required):
                raise ValueError("product phase is missing required arguments")
            values = {
                "alignment_mode": args.alignment_mode,
                "vulkan_backend": parse_bool(args.vulkan_backend),
                "compiled_kernels": parse_bool(args.compiled_kernels),
                "spirv_validation": parse_bool(args.spirv_validation),
                "native_viewer": parse_bool(args.native_viewer),
            }
            errors = product_errors(values)
        if errors:
            raise ValueError("; ".join(errors))
        print(f"{args.alignment_mode} {args.phase} build policy is satisfied")
        return 0
    except ValueError as error:
        parser.exit(1, f"build policy: {error}\n")


if __name__ == "__main__":
    sys.exit(main())
