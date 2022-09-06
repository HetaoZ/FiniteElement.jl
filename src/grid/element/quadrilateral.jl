

function init_volume(elem::Quadrilateral, nodes)
    x = elem_x0(elem, nodes)
    return polygon_area(x[1,:], x[2,:])
end

function volume(elem::Quadrilateral, nodes::Vector{Node})
    x = elem_x(elem, nodes)
    return polygon_area(x[1,:], x[2,:])
end

"""
Quad4's strain matrix
"""
function strain_matrix(::Quadrilateral, x, r, s)
    # x = elem_x(elem, nodes)
    x1 = x[1,1]; y1 = x[2,1]
    x2 = x[1,2]; y2 = x[2,2]
    x3 = x[1,3]; y3 = x[2,3]
    x4 = x[1,4]; y4 = x[2,4]

    Nr1, Nr2, Nr3, Nr4 = r-1, 1-r, r, -r
    Ns1, Ns2, Ns3, Ns4 = s-1, -s, s, 1-s

    J11 = Nr1*x1+Nr2*x2+Nr3*x3+Nr4*x4
    J12 = Nr1*y1+Nr2*y2+Nr3*y3+Nr4*y4
    J21 = Ns1*x1+Ns2*x2+Ns3*x3+Ns4*x4
    J22 = Ns1*y1+Ns2*y2+Ns3*y3+Ns4*y4

    detJ = J11*J22-J21*J12

    b11= J22*Nr1-J12*Ns1
    b21=-J21*Nr1+J11*Ns1

    b12= J22*Nr2-J12*Ns2
    b22=-J21*Nr2+J11*Ns2

    b13= J22*Nr3-J12*Ns3
    b23=-J21*Nr3+J11*Ns3

    b14= J22*Nr4-J12*Ns4
    b24=-J21*Nr4+J11*Ns4

    B = [b11  0    b12  0    b13  0    b14  0  ;
         0    b21  0    b22  0    b23  0    b24;
         b21  b11  b22  b12  b23  b13  b24  b14] 
    return B ./ detJ, detJ
end

function det_jacobi(::Quadrilateral, x)
    x1 = x[1,1]; y1 = x[2,1]
    x2 = x[1,2]; y2 = x[2,2]
    x3 = x[1,3]; y3 = x[2,3]
    x4 = x[1,4]; y4 = x[2,4]

    r, s = 0.5, 0.5
    Nr1, Nr2, Nr3, Nr4 = r-1, 1-r, r, -r
    Ns1, Ns2, Ns3, Ns4 = s-1, -s, s, 1-s

    J11 = Nr1*x1+Nr2*x2+Nr3*x3+Nr4*x4
    J12 = Nr1*y1+Nr2*y2+Nr3*y3+Nr4*y4
    J21 = Ns1*x1+Ns2*x2+Ns3*x3+Ns4*x4
    J22 = Ns1*y1+Ns2*y2+Ns3*y3+Ns4*y4

    detJ = J11*J22-J21*J12
    return detJ  
end