# Original retained control producer

These diagnostic tools build a private guarded-product candidate from preserved
complete shaders and run the eight original retained controls. Both the current
retained implementation and original audited test observer compile against the
same private generated header. All forty arrays are checked in the linked ELF.

Restore the files into `out/guarded-product-high-controls` inside the canonical
Sirius checkout and run there from the repository root. Use the existing
`linux-gcc` build. Preparation requires the local primary products and preserved
screen/observer evidence named in `prepare.py`; ignored binary evidence is not
backed up by this source tag. The checkout must remain clean at its stated
revision, and both diagnostic roots must be absent before preparation.

```sh
python3 -B out/guarded-product-high-controls/prepare.py
python3 -B out/guarded-product-high-controls/run_owned.py compile-compute-sealed
python3 -B out/guarded-product-high-controls/run_owned.py compile-observer
python3 -B out/guarded-product-high-controls/run_owned.py link-candidate
python3 -B out/guarded-product-high-controls/verify_producer.py
python3 -B out/guarded-product-high-controls/prepare_runtime.py
python3 -B out/guarded-product-high-controls/run_controls.py
python3 -B out/guarded-product-high-controls/check_controls.py
python3 -B out/guarded-product-high-controls/verify_control_routes.py
```

The controls retain the original fixtures, assertions, allocation cap and
300-second / 4-GiB RSS / 2-GiB reserve governors. The unchanged timestamp reader
checks complete readback, ordered marker and actual command/module records.
The separate Endpoint reader adds the original forced core route observations.
None of these commands grants performance acceptance, production adoption,
rendered-frame or full scientific/platform/release qualification.

Preservation and cleanup helpers apply to the recorded completed batch. They
require terminal owners, recoverable consumed input versions, exact archived
bytes/modes/topology, pushed source and fresh scoped consumer checks. Read their
preconditions before using them; do not execute them from tracked diagnostics.
