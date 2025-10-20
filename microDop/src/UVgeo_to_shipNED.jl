"""
Rotate U,V Earth-coordinate winds to ship North-East-Down coordinate
winds uNship (wind to forward), vEship (wind to starboard),
given the compass heading hed (degrees) of the ship. Compass heading 
North is 0, East is 90 degree.
""""
function UVgeo_to_shipNED(U,V, hed)
    # ship N is forward
    # ship E is to starboard
    #     (D is down)
    uNship =  V*cosd(hed) + U*sind(hed)
    vEship = -V*sind(hed) + U*cosd(hed)
    uNship, vEship
end
