#!/usr/bin/env python3
"""Separate post-owner byte reader; never launch a compiler, translator or provider."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import stat

ROOT = Path(__file__).absolute().parent


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def identity(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def frozen_inputs():
    freeze_path = ROOT / "input-freeze.json"
    freeze = json.loads(freeze_path.read_text())
    for name, pin in freeze["files"].items():
        path = ROOT / name
        require(not path.is_symlink() and path.is_file(), "missing/replaced frozen input: " + name)
        require({**identity(path), "mode": stat.S_IMODE(path.stat().st_mode)} == pin,
                "changed frozen input: " + name)
    for name, target in freeze["symlinks"].items():
        path = ROOT / name
        require(path.is_symlink() and os.readlink(path) == target and path.resolve().is_relative_to(ROOT),
                "changed/escaping frozen symlink: " + name)
    return identity(freeze_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stock", type=Path, required=True)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--result", type=Path, required=True)
    args = parser.parse_args()
    result_path = args.result.absolute()
    require(result_path.resolve().is_relative_to(ROOT / "observations") and not result_path.exists(),
            "result must be fresh inside owned observations")
    result = {"passed": False, "compiler_or_provider_launched": False, "qualification_claimed": False}
    try:
        result["input_freeze"] = frozen_inputs()
        expected_stock = (ROOT / "captured/initialize.metal").read_bytes()
        stock = args.stock.read_bytes()
        require(stock == expected_stock, "stock translation does not equal captured original MSL; reconstruction unresolved")
        result["stock"] = {"path": str(args.stock.absolute()), **identity(args.stock)}
        if args.candidate:
            candidate = args.candidate.read_bytes()
            expected = (ROOT / "captured/initialize-without-forced-inline.metal").read_bytes()
            require(candidate == expected, "candidate translation differs from the exact annotation-only control")
            require(candidate == stock.replace(b"static inline __attribute__((always_inline))", b"static inline"),
                    "candidate/stock relation is not annotation-only")
            result["candidate"] = {"path": str(args.candidate.absolute()), **identity(args.candidate)}
            result["decision"] = "TRANSLATION_BYTES_MATCH_ONLY_NATIVE_BUILD_PENDING"
        else:
            result["decision"] = "STOCK_TRANSLATION_MATCHES_CANDIDATE_TRANSLATION_PENDING"
        require(frozen_inputs() == result["input_freeze"], "freeze changed during separate reader")
        result["passed"] = True
    except Exception as error:
        result["error"] = f"{type(error).__name__}: {error}"
    result_path.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps(result))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
