## Initialize Librarys
rm(list = ls())
library(rgl)

## Return the Crossproduct of two vectors
#
getCrossProduct = function(x, y) {
    return(c(x[2]*y[3]-x[3]*y[2], -x[1]*y[3]+x[3]*y[1], x[1]*y[2]-x[2]*y[1]))
}

## Return the normalized Vector
#
normalize_vector <- function(x) {
    return(sqrt(sum(x^2)))
}

## Return the Rotationmatrix between to vector
#
getRotM <- function(x, y) {
    z = getCrossProduct(x, y) / normalize_vector(getCrossProduct(x, y))
    a = atan2(normalize_vector(getCrossProduct(x, y)), x %*% y)
    return(matrix(c(cos(a) + z[1] * z[1] * (1 - cos(a)), z[1] * z[2] * (1 - cos(a)) - z[3] * sin(a), z[1] * z[3] * (1 - cos(a)) + z[2] * sin(a),
                    z[1] * z[2] * (1 - cos(a)) + z[3] * sin(a),  cos(a) + z[2] * z[2] * (1 - cos(a)), z[2] * z[3] * (1 - cos(a)) - z[1] * sin(a),
                    z[1] * z[3] * (1 - cos(a)) - z[2] * sin(a),  z[2] * z[3] * (1 - cos(a)) + z[1] * sin(a), cos(a) + z[3] * z[3] * (1 - cos(a))), ncol=3, nrow=3, byrow = TRUE))
}

## Return the center of mass
##
centerofmass <- function(xyz) {
    com = rep(0, 3)
    for(n in 1:nrow(xyz)) {
        com[1] = com[1] + xyz[n, 1]
        com[2] = com[2] + xyz[n, 2]
        com[3] = com[3] + xyz[n, 3]
    }
    com[1] = com[1] / nrow(xyz)
    com[2] = com[2] / nrow(xyz)
    com[3] = com[3] / nrow(xyz)
    return(com);
}

## Return the model as a object from a xyz file
#
getModel <-function(url) {
    f=file(url)
    f_lines=strsplit(readLines(f),"\n")
    close(f)
    n=as.numeric(f_lines[[1]])
    m=matrix(nrow=n,ncol=3)
    l=c()
    for(i in 3:(n+2))
    {
        c=strsplit(f_lines[[i]],"[\t ]+")
        l=c(l,c[[1]][1])
        m[i-2,1]=as.numeric(c[[1]][2])
        m[i-2,2]=as.numeric(c[[1]][3])
        m[i-2,3]=as.numeric(c[[1]][4])
    }
    return(list(l=l, model=m));
}

## Return the Unitvector
#
getUnitVector <- function(x) {
    return(x/normalize_vector(x))
}

## Return the Radius of Gyration
#
radgyr <- function(xyz, com) {
    rsq = matrix(data = 0, nrow=nrow(xyz), ncol=ncol(xyz))
    for(n in 1:nrow(xyz)) {
        rsq[n, 1] = (xyz[n, 1] - com[1])^2
        rsq[n, 2] = (xyz[n, 2] - com[2])^2
        rsq[n, 3] = (xyz[n, 3] - com[3])^2
    }
    sqrs = rep(0, 4)
    sqrs[1] = sum(rsq)
    sqrs[2] = sum(rsq[, c(2,3)])
    sqrs[3] = sum(rsq[, c(1,3)])
    sqrs[4] = sum(rsq[, c(1,2)])
    rog_sq = sum(sqrs)/nrow(xyz)
    return(sqrt(rog_sq))
}

## Reset the center of mass plus entire obejct to a new location
#
resetcenter <- function(comstay, change) {
    comchange=centerofmass(change)
    xdist=abs(comstay[1]-comchange[1])
    ydist=abs(comstay[2]-comchange[2])
    zdist=abs(comstay[3]-comchange[3])
    if(comchange[1]>comstay[1]) { xdist=xdist-(2*xdist) }
    if(comchange[2]>comstay[2]) { ydist=ydist-(2*ydist) }
    if(comchange[3]>comstay[3]) { zdist=zdist-(2*zdist) }
    for(n in 1:nrow(change)) {
      change[n, 1] = change[n, 1] + xdist
      change[n, 2] = change[n, 2] + ydist
      change[n, 3] = change[n, 3] + zdist
    }
    return(change)
}

## Plot the entire alignment
#
plotAlignment <- function(m1, m2, l1, l2, iV1n, iV2n, iV3n, jV1n, jV2n, jV3n) {
    par3d("windowRect"= c(50,50,600,600))
    plot3d(m1,xlim=c(-20,20),ylim=c(-20,20),zlim=c(-20,20),box=F,axes=F,xlab="",ylab="",zlab="",col="red")
    lines3d(m1,col="red")
    text3d(m1,text=l1,col="red",cex=0.8)
    plot3d(m2,col="green",add=T)
    lines3d(m2,col="green")
    text3d(m2,text=l1,col="green",cex=0.8)
    lines3d(x=c(0, iV1n[1]),y=c(0, iV1n[2]),z=c(0, iV1n[3]), col="black")
    lines3d(x=c(0, iV2n[1]),y=c(0, iV2n[2]),z=c(0, iV2n[3]), col="black")
    lines3d(x=c(0, iV3n[1]),y=c(0, iV3n[2]),z=c(0, iV3n[3]), col="black")
    lines3d(x=c(0, jV1n[1]),y=c(0, jV1n[2]),z=c(0, jV1n[3]), col="purple")
    lines3d(x=c(0, jV2n[1]),y=c(0, jV2n[2]),z=c(0, jV2n[3]), col="purple")
    lines3d(x=c(0, jV3n[1]),y=c(0, jV3n[2]),z=c(0, jV3n[3]), col="purple")
    spheres3d(jV1n, radius = 0.05)
    spheres3d(jV2n, radius = 0.05)
    spheres3d(jV3n, radius = 0.05)
}

## Calculate the RMSD
#
rmsd <- function(v, w) {
    rmsd=0
    for(n in 1:nrow(v)) {
        rmsd=rmsd+((v[n, 1]-w[n, 1])^2+(v[n, 2]-w[n, 2])^2+(v[n, 3]-w[n, 3])^2)
    }
    return(sqrt(rmsd/nrow(v)))
}

## Rotate the entire model plus moment of inertia 
#
rotateEntireModel <- function(moi1, moi2, moi3, m, R) {
    moi1 <- as.vector(t(R %*% moi1))
    moi2 <- as.vector(t(R %*% moi2))
    moi3 <- as.vector(t(R %*% moi3))
    for(n in 1:nrow(m)) {
        x = c(m[n, 1], m[n, 2], m[n, 3])
        x <- R %*% x
        m[n, 1] = x[1]
        m[n, 2] = x[2]
        m[n, 3] = x[3]
    }
    return(list(moi1=moi1,moi2=moi2,moi3=moi3,model=m))
}

## export a model to a file
#
exportModel <- function (m, l, filename) {
    cat(nrow(m), file = filename, sep="\n")  # Apply cat & append
    for(n in 1:nrow(m)) {
        line= paste(l[n],m[n, 1],m[n, 2],m[n, 3],sep="\t")
        cat(line, file = filename, append = TRUE, sep="\n")
    }
}

## MAIN()
#
main <- function() {
    # Initialize Model Data
    file1 = getModel("struct/homework_02/1zta_ca_1.xyz")
    file2 = getModel("struct/homework_02/1zta_ca_2.xyz")
    m1 = file1$model
    m2 = file2$model
    l1 = file1$l
    l2 = file2$l
    
    # Initialize MOI Data -> sure here we can also directly read from file
    iV1=c(0.541059,-0.299595,0.785810)
    iV2=c(0.787317,-0.103191,-0.581439)
    iV3=c(0.010829,0.039588,0.007637)
    jV1=c(0.027761,-0.103066,0.994287)
    jV2=c(-0.010364,-0.974455,-0.100720)
    jV3=c(0.050197,-0.000385,-0.001441)
    iV1n=getUnitVector(iV1)
    iV2n=getUnitVector(iV2)
    iV3n=getUnitVector(iV3)
    jV1n=getUnitVector(jV1)
    jV2n=getUnitVector(jV2)
    jV3n=getUnitVector(jV3)

    # Actual center of mass
    com1 <- centerofmass(m1)
    com2 <- centerofmass(m2)
    print("Center of mass")
    print(com1)
    print(com2)

    # Radius of Gyration
    rgyr1 = radgyr(m1, com1)
    rgyr2 = radgyr(m2, com2)
    print("Radius of gyration")
    print(rgyr1)
    print(rgyr2)

    # Reset center of mass to c(0,0,0)
    m1new = resetcenter(c(0, 0, 0), m1)
    m2new = resetcenter(c(0, 0, 0), m2)
    com1new = centerofmass(m1new)
    com2new = centerofmass(m2new)

    # Rotate Model 2 to be superimposed to Model 1
    R1 = getRotM(jV1n, iV1n)
    R1sol=rotateEntireModel(jV1n, jV2n, jV3n, m2new, R1)
    R2 = getRotM(R1sol$moi2, iV2n)
    R2sol=rotateEntireModel(R1sol$moi1, R1sol$moi2, R1sol$moi3, R1sol$model, R2)
    R3 = getRotM(R2sol$moi3, iV3n)
    R3sol=rotateEntireModel(R2sol$moi1, R2sol$moi2, R2sol$moi3, R2sol$model, R3)
    jV1n <- R3sol$moi1
    jV2n <- R3sol$moi2
    jV3n <- R3sol$moi3
    m2new <- R3sol$model

    # Calculate RMSD
    print("RMSD")
    print(rmsd(m1new, m2new))

    # Plot final superimposed alignment
    plotAlignment(m1new, m2new, l1, l2, iV1n, iV2n, iV3n, jV1n, jV2n, jV3n)

    # export new models to xyz files
    exportModel(m1new, l1, "export_m1new.xyz")
    exportModel(m2new, l2, "export_m2new.xyz")
}

# start main function
main()