p=$1
cat <<HERE
{
    "dim": 2,
    "target": {
        "type": "orthotropic",
        "young":   [0.306026, 0.306026],
        "poisson": [0.154353, 0.154353],
        "shear":   [0.037147]
    },
    "target volume": 8.942765896271609,
    "initial_params": [0.8, $p],
    "radiusBounds": [0.1, 0.9],
    "translationBounds": [-0.15, 0.15],
    "blendingBounds": [1.1, 3.0]
}
HERE
