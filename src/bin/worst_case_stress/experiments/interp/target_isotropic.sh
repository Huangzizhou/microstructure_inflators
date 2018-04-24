young=$1
poisson=$2
cat <<END
{
    "dim": 2,
    "target": {
        "type": "isotropic",
        "young":   $young,
        "poisson": $poisson
    },
    "initial_params": [0.45, 0.5, 0.45, 0.45, 0.15, 0.15, 0.15, 0.15, 0.15, 0.01, 0.01, 0.01, 0.01, 0.01],
    "radiusBounds": [0.015, 0.2],
    "translationBounds": [0.3, 0.80],
    "blendingBounds": [0.0078125, 0.2]
}
END
