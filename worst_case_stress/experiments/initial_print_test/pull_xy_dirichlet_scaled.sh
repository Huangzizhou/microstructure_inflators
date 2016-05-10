magnitude=$1
cat <<END
{
    "no_rigid_motion": false,
    "regions": [
        { "type": "dirichlet", "value": [0, -$magnitude, 0], "box%": { "minCorner": [-0.0001, -0.0001, 0.0], "maxCorner": [1.0001, 0.0001, 0.0] } },
        { "type": "dirichlet", "value": [0,  $magnitude, 0], "box%": { "minCorner": [-0.0001,  0.9999, 0.0], "maxCorner": [1.0001, 1.0001, 0.0] } },
        { "type": "dirichlet", "value": [-$magnitude, 0, 0], "box%": { "minCorner": [-0.0001, -0.0001, 0.0], "maxCorner": [0.0001, 1.0001, 0.0] } },
        { "type": "dirichlet", "value": [ $magnitude, 0, 0], "box%": { "minCorner": [ 0.9999, -0.0001, 0.0], "maxCorner": [1.0001, 1.0001, 0.0] } }
    ]
}
END