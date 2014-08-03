// GMSH script to render the requested fields in a .msh file
// the angle bracket tags must be replaced appropriately before using (e.g. as
// done in make_flipper.pl)
General.SmallAxes = 0;
General.GraphicsWidth = 1024;
General.GraphicsHeight = 768;
Function DrawField
    For i In {0:<LAST_ITER>}
        name = StrCat(Sprintf("%g ", i ), field);
        For v In {0:PostProcessing.NbViews-1}
            If (StrCmp(View[v].Name, name) == 0)
                View[v].Visible = 1;
                // Vector field: displacement
                View[v].VectorType = 5;
                View[v].ShowElement = 1;
                View[v].Name = "";
                Draw;
                Print StrCat(StrCat(Sprintf("<PREFIX>.%g.", i), field),".png");
                View[v].Visible = 0;
            EndIf
        EndFor
    EndFor
Return

// Frustratingly, arrays of strings don't work in GMSH...
// Have PERL generate draw calls of the form:
// field = "u"; Call DrawField;
// Also, each of these has to be be on a separate line because of more GMSH
// fail.
<DRAW_CALLS>
Exit;
