

int main(int argc, char** argv) {
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " inpfile outfile  " << endl;
    return -1;
  }

  // open up the file to write to ...
  ifstream in(argv[1],std::ios::in|std::ios::binary);
  ofstream out(argv[2],std::ios::out|std::ios::binary);
  unsigned int n;
  in.read((char*)&n,sizeof(unsigned int));
  out.write((char*)&n,sizeof(unsigned int));

  double *x = new double[n];
  double *y = new double[n];
  double *z = new double[n];
  double maxX,minX,maxY,minY,maxZ,minZ;
  double meanX=0.0;
  double meanY=0.0;
  double meanZ=0.0;
  double inpMinX,inpMaxX,inpMinY,inpMaxY,inpMinZ,inpMaxZ;
  in.read((char*)&(x[0]),sizeof(double));
  in.read((char*)&(y[0]),sizeof(double));
  in.read((char*)&(z[0]),sizeof(double));
  inpMinX = inpMaxX = x[0];
  inpMinY = inpMaxY = y[0];
  inpMinZ = inpMaxZ = z[0];
  meanX += x[0];
  meanY += y[0];
  meanZ += z[0];
  for (unsigned int i = 1; i < n; i++) {
    in.read((char*)&(x[i]),sizeof(double));
    in.read((char*)&(y[i]),sizeof(double));
    in.read((char*)&(z[i]),sizeof(double));
    meanX += x[i];
    meanY += y[i];
    meanZ += z[i];
    if(inpMinX > x[i]) {
      inpMinX =  x[i];
    }
    if(inpMinY > y[i]) {
      inpMinY =  y[i];
    }
    if(inpMinZ > z[i]){
      inpMinZ =  z[i];
    }
    if(inpMaxX < x[i]) {
      inpMaxX = x[i];
    }
    if(inpMaxY < y[i]) {
      inpMaxY = y[i];
    }
    if(inpMaxZ < z[i]){
      inpMaxZ = z[i];
    }
  }
  meanX = meanX/((double)n);
  meanY = meanY/((double)n);
  meanZ = meanZ/((double)n);
  for(unsigned int i=0;i<n;i++) {
    x[i] = (x[i]-meanX)*2.0/(inpMaxX-inpMinX);
    y[i] = (y[i]-meanY)*2.0/(inpMaxY-inpMinY);
    z[i] = (z[i]-meanZ)*2.0/(inpMaxZ-inpMinZ);
  }

  double **xyz =  (double**)new char*[n];
  xyz[0] = new double[3];
  xyz[0][0] =  exp(x[0]);
  xyz[0][1] =  exp(y[0]);
  xyz[0][2] =  exp(z[0]);
  minX = maxX = xyz[0][0];
  minY = maxY = xyz[0][1];
  minZ = maxZ = xyz[0][2]; 
  for (unsigned int i = 1; i < n; i++) {
    xyz[i] = new double[3];
    xyz[i][0] =  exp(x[i]);
    xyz[i][1] =  exp(y[i]);
    xyz[i][2] =  exp(z[i]);
    if(minX > xyz[i][0]) {
      minX =  xyz[i][0];
    }
    if(minY > xyz[i][1]) {
      minY =  xyz[i][1];
    }
    if(minZ > xyz[i][2]){
      minZ =  xyz[i][2];
    }
    if(maxX < xyz[i][0]) {
      maxX = xyz[i][0];
    }
    if(maxY < xyz[i][1]) {
      maxY = xyz[i][1];
    }
    if(maxZ < xyz[i][2]){
      maxZ = xyz[i][2];
    }
  }
  for (unsigned int i = 0; i < n; i++) {
    xyz[i][0] = (xyz[i][0]-minX)/( (1.0 + 1.0e-3)*(maxX-minX) );
    xyz[i][1] = (xyz[i][1]-minY)/( (1.0 + 1.0e-3)*(maxY-minY) );
    xyz[i][2] = (xyz[i][2]-minZ)/( (1.0 + 1.0e-3)*(maxZ-minZ) );
    out.write((char*)xyz[i],3*sizeof(double));
    delete [] xyz[i];
  }
  delete [] xyz;
  delete [] x;
  delete [] y;
  delete [] z;
  in.close();
  out.close();
  return 0;
}

