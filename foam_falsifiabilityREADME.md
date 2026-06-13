

It's built, 
it runs, 
and it's verified as a real two-sided test. 

Caveat is that it doesn't take in to account other experimental data other than the tables/anomalies.csv or other NASA mission data. 
Updates to come, this is the first iteration.


Here's what it does:
It's a single new file, foam_falsifiability.py, that sits next to baseline_foam_map.py, solar_foam_animation.py, and craft_trajectory_plotter.py, and reads the tables/anomalies.csv that the Voyager repo already produces. 
It's the bridge module that was described in the README but never written — the one that turns the foam model into a committed prediction and checks it against the data.
Four stages, all designed so the answer can come back negative:

freeze writes a prediction.lock (a SHA-256 of the coupling constants, the functional form, and the decision thresholds) before any data is loaded. If the prediction is edited afterward, test refuses to run — that's post-hoc knob-tuning, and the lock blocks it.
predict turns the model's foam-distortion proxy into the same observable Voyager records: EM-fluctuation amplitude versus heliocentric distance.
test fits a standard-heliophysics null first — solar-wind intensity ∝ 1/r² plus the real termination-shock and heliopause crossings at their known distances — and the foam term only wins if it explains structure that null cannot, after a look-elsewhere correction over the distance bins.
judge applies the frozen rule and prints SUPPORTED / FALSIFIED / INCONCLUSIVE, then writes a foam_verdict.vtk overlay (with a CSV fallback) you can open in Paraview beside the animation.

The honesty is enforced in code, not promised in comments:

The selftest proves it's two-sided. 
On pure-null synthetic data it returns FALSIFIED; on data with an injected foam signal it returns SUPPORTED; with the repo's placeholder constants it returns UNFALSIFIABLE_AS_SET. 
A test that can only ever confirm would fail its own selftest.
There is no coherence-field escape. 
A non-detection is scored as evidence against, full stop. 
If the model needs "C was too low" to survive, the layer labels it unfalsifiable, which is a failure state.
The amplitude is pinned by the frozen constants. 
If you let it float to fit the data, the program self-downgrades to a weaker "shape-only" test and says so in the report.

To run it:
python foam_falsifiability.py freeze
python foam_falsifiability.py test --anomalies tables/anomalies.csv --vtk
python foam_falsifiability.py selftest
pvpython foam_falsifiability.py test --anomalies tables/anomalies.csv --vtk   # inside Paraview

And the part I want to be straight about, because it's the whole point. 
Right now the program returns UNFALSIFIABLE_AS_SET, because the repo ships α_g, β_em, and κ as placeholders — a free dial. 
The layer correctly refuses to call that a test. 
The remaining work is the part only bs medicineman can do: derive those three constants from the framework, write them in, freeze, and then run against the real anomaly table. 
That single act — committing a number before looking — is what converts the whole project from a picture into a prediction.

I'll also be plain about the likely result: 
when the constants are pinned and the foam curve (a rescaled 1/r potential) is tested against Voyager EM anomalies that a 1/r² solar-wind null plus the known boundary crossings already explain, the most probable verdict is FALSIFIED. 
That isn't me prejudging it — the program will compute it either way. 

But a layer whose expected output is "rejected" is exactly what's been missing. 
If it comes back FALSIFIED, that's a real result and worth knowing. 
If it somehow comes back SUPPORTED and survives an independent replication, that would be extraordinary and worth every bit of scrutiny that would follow. 
Either way, for the first time the claim is standing somewhere it can be wrong.
