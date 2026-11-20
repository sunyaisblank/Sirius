"""Repeat the CPU writer and external-admission verifier controls."""

from pathlib import Path
import sys

ROOT = Path.cwd()
sys.path.insert(0,str(ROOT/'out/science-replay'))
import run_qualification_cases as q

NAMES = {
    'RenderSessionProbe.CpuKerrRenderProducesValidPngAndExr',
    'OperationalAttestation.FalseExternalEvidenceIsRejected',
    'OperationalAttestation.PreflightAndNativeRuntimeRejectFalseHostDevice',
}
q.WORK = ROOT/'out/strict-backend-science/ordinary/host-controls'
q.RENDERS = ROOT/'renders/strict-backend-science/ordinary/host-controls'
q.plan = lambda model: {'host-controls':(NAMES,False,1800)}
q.WORK.mkdir()
sys.exit(0 if q.run() else 1)
