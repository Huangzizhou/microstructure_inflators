blend=$1
cat <<HERE
{
    "dim": 2,
    "target": {
        "type": "orthotropic",
        "young":   [2.933545, 2.933545],
        "poisson": [0.27186, 0.27186],
        "shear":   [0.87212]
    },
    "initial_params": [0.75, 0.75, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, $blend, $blend, $blend, $blend, $blend],
    "radiusBounds": [0.05, 0.2],
    "translationBounds": [0.25, 1.0],
    "blendingBounds": [8, 128]
}
HERE
