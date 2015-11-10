//#include <stdint.h>
//#include "dendro_vtk.h"
//
//const int CHARS_PER_LINE = 72;
//
//#ifndef DENDRO_VTK_DOUBLES
//#define DENDRO_VTK_FLOAT_NAME "Float32"
//#define DENDRO_VTK_FLOAT_TYPE float
//#else
//#define DENDRO_VTK_FLOAT_NAME "Float64"
//#define DENDRO_VTK_FLOAT_TYPE double
//#endif
//
//#ifndef DENDRO_VTK_BINARY
//#define DENDRO_VTK_ASCII 1
//#define DENDRO_VTK_FORMAT_STRING "ascii"
//#else
//#define DENDRO_VTK_FORMAT_STRING "binary"
//#endif
//
//void base64_init_encodestate(base64_encodestate* state_in)
//{
//  state_in->step = step_A;
//  state_in->result = 0;
//  state_in->stepcount = 0;
//}
//
//char base64_encode_value(char value_in)
//{
//  static const char* encoding = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
//  if (value_in > 63) return '=';
//  return encoding[(int)value_in];
//}
//
//int base64_encode_block(const char* plaintext_in, int length_in, char* code_out, base64_encodestate* state_in)
//{
//  const char* plainchar = plaintext_in;
//  const char* const plaintextend = plaintext_in + length_in;
//  char* codechar = code_out;
//  char result;
//  char fragment;
//
//  result = state_in->result;
//
//  switch (state_in->step)
//  {
//    while (1)
//    {
//      case step_A:
//        if (plainchar == plaintextend)
//        {
//          state_in->result = result;
//          state_in->step = step_A;
//          return codechar - code_out;
//        }
//      fragment = *plainchar++;
//      result = (fragment & 0x0fc) >> 2;
//      *codechar++ = base64_encode_value(result);
//      result = (fragment & 0x003) << 4;
//      case step_B:
//        if (plainchar == plaintextend)
//        {
//          state_in->result = result;
//          state_in->step = step_B;
//          return codechar - code_out;
//        }
//      fragment = *plainchar++;
//      result |= (fragment & 0x0f0) >> 4;
//      *codechar++ = base64_encode_value(result);
//      result = (fragment & 0x00f) << 2;
//      case step_C:
//        if (plainchar == plaintextend)
//        {
//          state_in->result = result;
//          state_in->step = step_C;
//          return codechar - code_out;
//        }
//      fragment = *plainchar++;
//      result |= (fragment & 0x0c0) >> 6;
//      *codechar++ = base64_encode_value(result);
//      result  = (fragment & 0x03f) >> 0;
//      *codechar++ = base64_encode_value(result);
//
//      ++(state_in->stepcount);
//      if (state_in->stepcount == CHARS_PER_LINE/4)
//      {
//        *codechar++ = '\n';
//        state_in->stepcount = 0;
//      }
//    }
//  }
//  /* control should not reach here */
//  return codechar - code_out;
//}
//
//int base64_encode_blockend(char* code_out, base64_encodestate* state_in)
//{
//  char* codechar = code_out;
//
//  switch (state_in->step)
//  {
//    case step_B:
//      *codechar++ = base64_encode_value(state_in->result);
//      *codechar++ = '=';
//      *codechar++ = '=';
//      break;
//    case step_C:
//      *codechar++ = base64_encode_value(state_in->result);
//      *codechar++ = '=';
//      break;
//    case step_A:
//      break;
//  }
//  *codechar++ = '\n';
//
//  return codechar - code_out;
//}
//
//
//int
//sc_vtk_write_binary (FILE * vtkfile, char *numeric_data, size_t byte_length)
//{
//  size_t              chunks, chunksize, remaining, writenow;
//  size_t              code_length, base_length;
//  unsigned int        int_header;
//  char               *base_data;
//  base64_encodestate  encode_state;
//
//  /* VTK format used 32bit header info */
//  assert (byte_length <= (size_t) UINT32_MAX);
//
//  /* This value may be changed although this is not tested with VTK */
//  chunksize = (size_t) 1 << 15; /* 32768 */
//  int_header = (uint32_t) byte_length;
//
//  /* Allocate sufficient memory for base64 encoder */
//  code_length = 2 * std::max (chunksize, sizeof (int_header));
//  code_length = std::max (code_length, 4) + 1;
//  base_data = new char[code_length];
//
//  base64_init_encodestate (&encode_state);
//  base_length = base64_encode_block ((char *) &int_header, sizeof (int_header), base_data, &encode_state);
//  assert (base_length < code_length);
//  base_data[base_length] = '\0';
//  (void) fwrite (base_data, 1, base_length, vtkfile);
//
//  chunks = 0;
//  remaining = byte_length;
//  while (remaining > 0) {
//    writenow = std::min (remaining, chunksize);
//    base_length = base64_encode_block (numeric_data + chunks * chunksize,
//                                       writenow, base_data, &encode_state);
//    assert (base_length < code_length);
//    base_data[base_length] = '\0';
//    (void) fwrite (base_data, 1, base_length, vtkfile);
//    remaining -= writenow;
//    ++chunks;
//  }
//
//  base_length = base64_encode_blockend (base_data, &encode_state);
//  assert (base_length < code_length);
//  base_data[base_length] = '\0';
//  (void) fwrite (base_data, 1, base_length, vtkfile);
//
//  delete [] base_data;
//  if (ferror (vtkfile)) {
//    return -1;
//  }
//  return 0;
//}
//
//void
//dendro_vtk_write_file (std::vector<ot::TreeNode>& tree, const char *baseName)
//{
//  int                 retval;
//  int                 dendro_vtk_write_rank = 1;
//
//  retval = dendro_vtk_write_header (tree, baseName, NULL, NULL, dendro_vtk_write_rank);
//  if (!retval) throw std::ios_base::failure("VTK: write header");
//  retval = dendro_vtk_write_footer (tree, baseName);
//  if (!retval) throw std::ios_base::failure("VTK: write footer");
//}
//
//int
//dendro_vtk_write_header (std::vector<ot::TreeNode> &tree, const char *baseName,
//                         const char *pointscalars, const char *pointvectors,
//                         int write_rank)
//{
//  return dendro_vtk_write_header_all (tree, baseName, pointscalars,
//                                      pointvectors, write_rank, 0, 0);
//}
//
//int
//dendro_vtk_write_header_wrap (std::vector<ot::TreeNode> &tree, const char *baseName,
//                              const char *pointscalars,
//                              const char *pointvectors, int write_rank,
//                              int rank_wrap)
//{
//  return dendro_vtk_write_header_all (tree, baseName, pointscalars,
//                                      pointvectors, write_rank, rank_wrap, 0);
//}
//
//int
//dendro_vtk_write_header_all (std::vector<ot::TreeNode> &tree, const char *baseName,
//                             const char *pointscalars,
//                             const char *pointvectors, int write_rank,
//                             int rank_wrap, int write_tag)
//{
//  double             *Xd = ma->X->e[0];
//  double             *Yd = ma->Y->e[0];
//  double             *Zd = ma->Z->e[0];
//  const int           mpirank = ma->mpirank;
//  const int           dim = ma->D;
//  const int           N = ma->N;
//  const int           Np = ma->Np;
//  const int           Nrp = ma->Nrp;
//  const int           Nvertices = ma->Nvertices;
//  const int           N3 = N * N * N;
//  const mangll_locidx_t K = ma->mesh->K;
//  const mangll_locidx_t Ncells = K * N3;
//  const mangll_locidx_t Ntotal = K * Np;
//
//  #ifdef MANGLL_VTK_ASCII
//  mangll_locidx_t     sk;
//#else
//  int                 retval;
//  uint8_t            *uint8_data;
//  mangll_tag_t       *tagt_data;
//  DendroIntL         *locidx_data;
//  DENDRO_VTK_FLOAT_TYPE *float_data;
//#endif
//
//  int                 i, j, k;
//  int                 mpiwrap;
//  mangll_locidx_t     ik, il;
//  char                vtufilename[BUFSIZ];
//  FILE               *vtufile;
//
//  /* Have each proc write to its own file */
//  snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", baseName, mpirank);
//  vtufile = fopen (vtufilename, "w");
//
//  if (vtufile == NULL) {
//    printf("Could not open %s for output!\n", vtufilename);
//    return -1;
//  }
//
//  fprintf (vtufile, "<?xml version=\"1.0\"?>\n");
//  fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
//#if defined(MANGLL_VTK_BINARY) && defined(MANGLL_VTK_COMPRESSION)
//  fprintf (vtufile, " compressor=\"vtkZLibDataCompressor\"");
//#endif
//#ifdef SC_WORDS_BIGENDIAN
//  fprintf (vtufile, " byte_order=\"BigEndian\">\n");
//#else
//  fprintf (vtufile, " byte_order=\"LittleEndian\">\n");
//#endif
//  fprintf (vtufile, "  <UnstructuredGrid>\n");
//  fprintf (vtufile,
//           "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
//           (long long) Ntotal, (long long) Ncells);
//  fprintf (vtufile, "      <Points>\n");
//
//  /* write point position data */
//  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
//               " NumberOfComponents=\"3\" format=\"%s\">\n",
//           DENDRO_VTK_FLOAT_NAME, DENDRO_VTK_FORMAT_STRING);
//
//#ifdef DENDRO_VTK_ASCII
//  for (il = 0; il < Ntotal; ++il) {
//#ifdef DENDRO_VTK_DOUBLES
//    fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", Xd[il], Yd[il], Zd[il]);
//#else
//    fprintf (vtufile, "          %16.8e %16.8e %16.8e\n", Xd[il], Yd[il], Zd[il]);
//#endif
//  }
//#else
//  float_data = new DENDRO_VTK_FLOAT_TYPE[dim * Ntotal];
//  for (il = 0; il < Ntotal; ++il) {
//    float_data[3 * il + 0] = (DENDRO_VTK_FLOAT_TYPE) Xd[il];
//    float_data[3 * il + 1] = (DENDRO_VTK_FLOAT_TYPE) Yd[il];
//    float_data[3 * il + 2] = (DENDRO_VTK_FLOAT_TYPE) Zd[il];
//  }
//
//  fprintf (vtufile, "          ");
//  retval = dendro_vtk_write_binary (vtufile, (char *) float_data, sizeof (*float_data) * dim * Ntotal);
//  fprintf (vtufile, "\n");
//  if (retval) {
//    printf("mangll_vtk: Error encoding points\n");
//    fclose (vtufile);
//    return -1;
//  }
//  delete [] float_data;
//#endif
//  fprintf (vtufile, "        </DataArray>\n");
//  fprintf (vtufile, "      </Points>\n");
//  fprintf (vtufile, "      <Cells>\n");
//
//  /* write connectivity data */
//  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"connectivity\""
//      " format=\"%s\">\n", DENDRO_VTK_LOCIDX, DENDRO_VTK_FORMAT_STRING);
//  assert (Nvertices == 8);
//#ifdef MANGLL_VTK_ASCII
//  for (ik = 0; ik < K; ++ik) {
//    for (k = 0; k < N; ++k) {
//      for (j = 0; j < N; ++j) {
//        for (i = 0; i < N; ++i) {
//          fprintf (vtufile,
//                   "          %lld %lld %lld %lld %lld %lld %lld %lld\n",
//                   (long long) ik * Np + k * Nrp * Nrp + j * Nrp + i,
//                   (long long) ik * Np + k * Nrp * Nrp + j * Nrp + (i + 1),
//                   (long long) ik * Np + k * Nrp * Nrp + (j + 1) * Nrp + (i +
//                                                                          1),
//                   (long long) ik * Np + k * Nrp * Nrp + (j + 1) * Nrp + i,
//                   (long long) ik * Np + (k + 1) * Nrp * Nrp + j * Nrp + i,
//                   (long long) ik * Np + (k + 1) * Nrp * Nrp + j * Nrp + (i +
//                                                                          1),
//                   (long long) ik * Np + (k + 1) * Nrp * Nrp + (j + 1) * Nrp +
//                   (i + 1),
//                   (long long) ik * Np + (k + 1) * Nrp * Nrp + (j + 1) * Nrp +
//                   i);
//        }
//      }
//    }
//  }
//#else
//  locidx_data = new DendroIntL[Ncells * Nvertices];
//  for (ik = 0, il = 0; ik < K; ++ik) {
//    for (k = 0; k < N; ++k) {
//      for (j = 0; j < N; ++j) {
//        for (i = 0; i < N; ++i) {
//          locidx_data[il + 0] = ik * Np + k * Nrp * Nrp + j * Nrp + i;
//          locidx_data[il + 1] = ik * Np + k * Nrp * Nrp + j * Nrp + (i + 1);
//          locidx_data[il + 2] =
//              ik * Np + k * Nrp * Nrp + (j + 1) * Nrp + (i + 1);
//          locidx_data[il + 3] = ik * Np + k * Nrp * Nrp + (j + 1) * Nrp + i;
//          locidx_data[il + 4] = ik * Np + (k + 1) * Nrp * Nrp + j * Nrp + i;
//          locidx_data[il + 5] = ik * Np + (k + 1) * Nrp * Nrp + j * Nrp + (i +
//                                                                           1);
//          locidx_data[il + 6] =
//              ik * Np + (k + 1) * Nrp * Nrp + (j + 1) * Nrp + (i + 1);
//          locidx_data[il + 7] =
//              ik * Np + (k + 1) * Nrp * Nrp + (j + 1) * Nrp + i;
//          il += 8;
//        }
//      }
//    }
//  }
//
//  fprintf (vtufile, "          ");
//  retval = dendro_vtk_write_binary (vtufile, (char *) locidx_data,
//                                    sizeof (*locidx_data) * Ncells *
//                                    Nvertices);
//  fprintf (vtufile, "\n");
//  if (retval) {
//    printf ("mangll_vtk: Error encoding connectivity\n");
//    fclose (vtufile);
//    return -1;
//  }
//#endif
//  fprintf (vtufile, "        </DataArray>\n");
//
//  /* write offset data */
//  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"offsets\""
//      " format=\"%s\">\n", MANGLL_VTK_LOCIDX, DENDRO_VTK_FORMAT_STRING);
//#ifdef DENDRO_VTK_ASCII
//  fprintf (vtufile, "         ");
//  for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
//    fprintf (vtufile, " %lld", (long long) (Nvertices * il));
//    if (!(sk % 8) && il != Ncells)
//      fprintf (vtufile, "\n         ");
//  }
//  fprintf (vtufile, "\n");
//#else
//  for (il = 1; il <= Ncells; ++il)
//    locidx_data[il - 1] = Nvertices * il;
//
//  fprintf (vtufile, "          ");
//  retval = dendro_vtk_write_binary (vtufile, (char *) locidx_data, sizeof (*locidx_data) * Ncells);
//  fprintf (vtufile, "\n");
//  if (retval) {
//    printf ("mangll_vtk: Error encoding offsets\n");
//    fclose (vtufile);
//    return -1;
//  }
//#endif
//  fprintf (vtufile, "        </DataArray>\n");
//
//  /* write type data */
//  fprintf (vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
//      " format=\"%s\">\n", DENDRO_VTK_FORMAT_STRING);
//#ifdef DENDRO_VTK_ASCII
//  fprintf (vtufile, "         ");
//  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
//    fprintf (vtufile, " 12");   /* 12:VTK_HEXAHEDRON */
//    if (!(sk % 20) && il != (Ncells - 1))
//      fprintf (vtufile, "\n         ");
//  }
//  fprintf (vtufile, "\n");
//#else
//  uint8_data = new uint8_t[Ncells];
//  for (il = 0; il < Ncells; ++il)
//    uint8_data[il] = 12;
//
//  fprintf (vtufile, "          ");
//  retval = dendro_vtk_write_binary (vtufile, (char *) uint8_data, sizeof (*uint8_data) * Ncells);
//  delete [] uint8_data;
//
//  fprintf (vtufile, "\n");
//  if (retval) {
//    printf ("dendro_vtk: Error encoding types\n");
//    fclose (vtufile);
//    return -1;
//  }
//#endif
//  fprintf (vtufile, "        </DataArray>\n");
//  fprintf (vtufile, "      </Cells>\n");
//
//  if (write_rank || rank_wrap > 0 || write_tag) {
//    char                linebuf[BUFSIZ];
//
//    i = 0;
//    j = snprintf (linebuf, BUFSIZ, "      <CellData Scalars=\"");
//    if (write_rank) {
//      j += snprintf (linebuf + j, BUFSIZ - j, "%smpirank", i++ ? "," : "");
//    }
//    if (rank_wrap > 0) {
//      j += snprintf (linebuf + j, BUFSIZ - j, "%smpiwrap", i++ ? "," : "");
//    }
//    if (write_tag) {
//      j += snprintf (linebuf + j, BUFSIZ - j, "%setag", i++ ? "," : "");
//    }
//    assert (i > 0);
//    fprintf (vtufile, "%s%s\n", linebuf, "\">");
//  }
//  if (write_rank) {
//    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"mpirank\""
//                 " format=\"%s\">\n", MANGLL_VTK_LOCIDX,
//             DENDRO_VTK_FORMAT_STRING);
//#ifdef DENDRO_VTK_ASCII
//    fprintf (vtufile, "         ");
//    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
//      fprintf (vtufile, " %d", mpirank);
//      if (!(sk % 20) && il != (Ncells - 1))
//        fprintf (vtufile, "\n         ");
//    }
//    fprintf (vtufile, "\n");
//#else
//    for (il = 0; il < Ncells; ++il)
//      locidx_data[il] = (mangll_locidx_t) mpirank;
//
//    fprintf (vtufile, "          ");
//    retval = dendro_vtk_write_binary (vtufile, (char *) locidx_data,
//                                      sizeof (*locidx_data) * Ncells);
//    fprintf (vtufile, "\n");
//    if (retval) {
//      printf ("dendro_vtk: Error encoding mpirank\n");
//      fclose (vtufile);
//      return -1;
//    }
//#endif
//    fprintf (vtufile, "        </DataArray>\n");
//  }
//  if (rank_wrap > 0) {
//    assert (rank_wrap <= ma->mpisize);
//    mpiwrap = (mpirank * rank_wrap) % ma->mpisize;
//    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"mpiwrap\""
//                 " format=\"%s\">\n", MANGLL_VTK_LOCIDX,
//             DENDRO_VTK_FORMAT_STRING);
//#ifdef DENDRO_VTK_ASCII
//    fprintf (vtufile, "         ");
//    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
//      fprintf (vtufile, " %d", mpiwrap);
//      if (!(sk % 20) && il != (Ncells - 1))
//        fprintf (vtufile, "\n         ");
//    }
//    fprintf (vtufile, "\n");
//#else
//    for (il = 0; il < Ncells; ++il)
//      locidx_data[il] = (mangll_locidx_t) mpiwrap;
//
//    fprintf (vtufile, "          ");
//    retval = dendro_vtk_write_binary (vtufile, (char *) locidx_data,
//                                      sizeof (*locidx_data) * Ncells);
//    fprintf (vtufile, "\n");
//    if (retval) {
//      printf ("dendro_vtk: Error encoding mpiwrap\n");
//      fclose (vtufile);
//      return -1;
//    }
//#endif
//    fprintf (vtufile, "        </DataArray>\n");
//  }
//  if (write_tag) {
//    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"etag\""
//        " format=\"%s\">\n", MANGLL_VTK_TAGT, DENDRO_VTK_FORMAT_STRING);
//#ifdef DENDRO_VTK_ASCII
//    fprintf (vtufile, "         ");
//    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
//      fprintf (vtufile, " %ld", (long) ma->mesh->EtoT[il / N3]);
//      if (!(sk % 20) && il != (Ncells - 1))
//        fprintf (vtufile, "\n         ");
//    }
//    fprintf (vtufile, "\n");
//#else
//    tagt_data = MANGLL_ALLOC (mangll_tag_t, Ncells);
//    for (il = 0; il < Ncells; ++il) {
//      tagt_data[il] = ma->mesh->EtoT[il / N3];
//    }
//    fprintf (vtufile, "          ");
//    retval = mangll_vtk_write_binary (vtufile, (char *) tagt_data,
//                                      sizeof (*tagt_data) * Ncells);
//    fprintf (vtufile, "\n");
//    if (retval) {
//      printf("dendro_vtk: Error encoding element tags\n");
//      fclose (vtufile);
//      return -1;
//    }
//    MANGLL_FREE (tagt_data);
//#endif
//    fprintf (vtufile, "        </DataArray>\n");
//  }
//  if (write_rank || rank_wrap > 0 || write_tag) {
//    fprintf (vtufile, "      </CellData>\n");
//  }
//#ifndef DENDRO_VTK_ASCII
//  delete  [] locidx_data;
//#endif
//
//  fprintf (vtufile, "      <PointData");
//  if (pointscalars != NULL)
//    fprintf (vtufile, " Scalars=\"%s\"", pointscalars);
//  if (pointvectors != NULL)
//    fprintf (vtufile, " Vectors=\"%s\"", pointvectors);
//  fprintf (vtufile, ">\n");
//
//  if (ferror (vtufile)) {
//    printf ("dendro_vtk: Error writing header\n");
//    fclose (vtufile);
//    return -1;
//  }
//  if (fclose (vtufile)) {
//    printf ("dendro_vtk: Error closing header\n");
//    return -1;
//  }
//  vtufile = NULL;
//
//  /* Only have the root write to the parallel vtk file */
//  if (mpirank == 0) {
//    char                pvtufilename[BUFSIZ];
//    FILE               *pvtufile;
//    snprintf (pvtufilename, BUFSIZ, "%s.pvtu", baseName);
//
//    pvtufile = fopen (pvtufilename, "w");
//    if (!pvtufile) {
//      printf ("Could not open %s for output!\n", vtufilename);
//      return -1;
//    }
//
//    fprintf (pvtufile, "<?xml version=\"1.0\"?>\n");
//    fprintf (pvtufile, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
//#if defined(DENDRO_VTK_BINARY) && defined(DENDRO_VTK_COMPRESSION)
//    fprintf (pvtufile, " compressor=\"vtkZLibDataCompressor\"");
//#endif
//#ifdef DENDRO_WORDS_BIGENDIAN
//    fprintf (pvtufile, " byte_order=\"BigEndian\">\n");
//#else
//    fprintf (pvtufile, " byte_order=\"LittleEndian\">\n");
//#endif
//
//    fprintf (pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
//    fprintf (pvtufile, "    <PPoints>\n");
//    fprintf (pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\""
//                 " NumberOfComponents=\"3\" format=\"%s\"/>\n",
//             DENDRO_VTK_FLOAT_NAME, DENDRO_VTK_FORMAT_STRING);
//    fprintf (pvtufile, "    </PPoints>\n");
//
//    if (write_rank || rank_wrap > 0 || write_tag) {
//      char                linebuf[BUFSIZ];
//
//      i = 0;
//      j = snprintf (linebuf, BUFSIZ, "      <PCellData Scalars=\"");
//      if (write_rank) {
//        j += snprintf (linebuf + j, BUFSIZ - j, "%smpirank", i++ ? "," : "");
//      }
//      if (rank_wrap > 0) {
//        j += snprintf (linebuf + j, BUFSIZ - j, "%smpiwrap", i++ ? "," : "");
//      }
//      if (write_tag) {
//        j += snprintf (linebuf + j, BUFSIZ - j, "%setag", i++ ? "," : "");
//      }
//      assert (i > 0);
//      fprintf (pvtufile, "%s%s\n", linebuf, "\">");
//    }
//    if (write_rank) {
//      fprintf (pvtufile,
//               "      <PDataArray type=\"%s\" Name=\"mpirank\" format"
//                   "=\"%s\"/>\n", MANGLL_VTK_LOCIDX, DENDRO_VTK_FORMAT_STRING);
//    }
//    if (rank_wrap > 0) {
//      fprintf (pvtufile,
//               "      <PDataArray type=\"%s\" Name=\"mpiwrap\" format="
//                   "\"%s\"/>\n", MANGLL_VTK_LOCIDX, DENDRO_VTK_FORMAT_STRING);
//    }
//    if (write_tag) {
//      fprintf (pvtufile,
//               "      <PDataArray type=\"%s\" Name=\"etag\" format="
//                   "\"%s\"/>\n", MANGLL_VTK_TAGT, DENDRO_VTK_FORMAT_STRING);
//    }
//    if (write_rank || rank_wrap > 0 || write_tag) {
//      fprintf (pvtufile, "    </PCellData>\n");
//    }
//
//    fprintf (pvtufile, "    <PPointData");
//    if (pointscalars != NULL)
//      fprintf (pvtufile, " Scalars=\"%s\"", pointscalars);
//    if (pointvectors != NULL)
//      fprintf (pvtufile, " Vectors=\"%s\"", pointvectors);
//    fprintf (pvtufile, ">\n");
//
//    if (ferror (pvtufile)) {
//      printf ("dendro_vtk: Error writing parallel header\n");
//      fclose (pvtufile);
//      return -1;
//    }
//    if (fclose (pvtufile)) {
//      printf ("dendro_vtk: Error closing parallel header\n");
//      return -1;
//    }
//  }
//
//  return 0;
//}
//
//int
//mangll_vtk_write_point_scalar (mangll_t * ma,
//                               const char *baseName,
//                               double *scalard, int sskip,
//                               const char *scalarName)
//{
//  char                vtufilename[BUFSIZ];
//  int                 procRank = ma->mpirank;
//  const int           Np = ma->Np;
//  const mangll_locidx_t K = ma->mesh->K;
//  const mangll_locidx_t Ntotal = K * Np;
//  mangll_locidx_t     il;
//  FILE               *vtufile;
//#ifndef MANGLL_VTK_ASCII
//  MANGLL_VTK_FLOAT_TYPE *float_data;
//  int                 retval;
//#endif
//
//  /* Have each proc write to its own file */
//  snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", baseName, procRank);
//  /* To be able to fseek in a file you can not open in append mode.
//   * so you need to open with r+ and fseek to SEEK_END.
//   */
//  vtufile = fopen (vtufilename, "r+");
//  if (vtufile == NULL) {
//    MANGLL_LERRORF ("Could not open %s for output!\n", vtufilename);
//    return -1;
//  }
//  fseek (vtufile, 0L, SEEK_END);
//
//  /* write point position data */
//  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"%s\""
//               " format=\"%s\">\n",
//           MANGLL_VTK_FLOAT_NAME, scalarName, MANGLL_VTK_FORMAT_STRING);
//
//#ifdef MANGLL_VTK_ASCII
//  for (il = 0; il < Ntotal; ++il) {
//#ifdef MANGLL_VTK_DOUBLES
//    fprintf (vtufile, "     %24.16e\n", scalard[sskip * il]);
//#else
//    fprintf (vtufile, "          %16.8e\n", scalard[sskip * il]);
//#endif
//  }
//#else
//  /* TODO Don't allocate the full size of the array, only allocate
//   * the chunk that will be passed to zlib and do this a chunk
//   * at a time.
//   */
//  float_data = MANGLL_ALLOC (MANGLL_VTK_FLOAT_TYPE, Ntotal);
//  for (il = 0; il < Ntotal; ++il) {
//    float_data[il] = (MANGLL_VTK_FLOAT_TYPE) scalard[sskip * il];
//  }
//
//  fprintf (vtufile, "          ");
//  retval = mangll_vtk_write_binary (vtufile, (char *) float_data,
//                                    sizeof (*float_data) * Ntotal);
//  fprintf (vtufile, "\n");
//  if (retval) {
//    MANGLL_LERROR ("mangll_vtk: Error encoding points\n");
//    fclose (vtufile);
//    return -1;
//  }
//  MANGLL_FREE (float_data);
//#endif
//  fprintf (vtufile, "        </DataArray>\n");
//
//  if (ferror (vtufile)) {
//    MANGLL_LERROR ("mangll_vtk: Error writing point scalar\n");
//    fclose (vtufile);
//    return -1;
//  }
//  if (fclose (vtufile)) {
//    MANGLL_LERROR ("mangll_vtk: Error closing point scalar\n");
//    return -1;
//  }
//  vtufile = NULL;
//
//  /* Only have the root write to the parallel vtk file */
//  if (procRank == 0) {
//    char                pvtufilename[BUFSIZ];
//    FILE               *pvtufile;
//    snprintf (pvtufilename, BUFSIZ, "%s.pvtu", baseName);
//
//    pvtufile = fopen (pvtufilename, "a");
//    if (!pvtufile) {
//      MANGLL_LERRORF ("Could not open %s for output!\n", vtufilename);
//      return -1;
//    }
//
//    fprintf (pvtufile, "      <PDataArray type=\"%s\" Name=\"%s\""
//                 " format=\"%s\"/>\n",
//             MANGLL_VTK_FLOAT_NAME, scalarName, MANGLL_VTK_FORMAT_STRING);
//
//    if (ferror (pvtufile)) {
//      MANGLL_LERROR ("mangll_vtk: Error writing parallel point scalar\n");
//      fclose (pvtufile);
//      return -1;
//    }
//    if (fclose (pvtufile)) {
//      MANGLL_LERROR ("mangll_vtk: Error closing parallel point scalar\n");
//      return -1;
//    }
//  }
//
//  return 0;
//}
//
//int
//mangll_vtk_write_point_vector (mangll_t * ma, const char *baseName,
//                               double *vd, int vskip, const char *vectorName)
//{
//  mangll_vtk_write_point_vector_123 (ma, baseName, vd + 0, vd + 1, vd + 2,
//                                     vskip, vectorName);
//
//  return 0;
//}
//
//int
//mangll_vtk_write_point_vector_123 (mangll_t * ma, const char *baseName,
//                                   double *v1d, double *v2d, double *v3d,
//                                   int vskip, const char *vectorName)
//{
//  char                vtufilename[BUFSIZ];
//  int                 procRank = ma->mpirank;
//  const int           dim = ma->D;
//  const int           Np = ma->Np;
//  const mangll_locidx_t K = ma->mesh->K;
//  const mangll_locidx_t Ntotal = K * Np;
//  mangll_locidx_t     il;
//  FILE               *vtufile;
//#ifndef MANGLL_VTK_ASCII
//  MANGLL_VTK_FLOAT_TYPE *float_data;
//  int                 retval;
//#endif
//
//  /* Have each proc write to its own file */
//  snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", baseName, procRank);
//  /* To be able to fseek in a file you can not open in append mode.
//   * so you need to open with r+ and fseek to SEEK_END.
//   */
//  vtufile = fopen (vtufilename, "r+");
//  if (vtufile == NULL) {
//    MANGLL_LERRORF ("Could not open %s for output!\n", vtufilename);
//    return -1;
//  }
//  fseek (vtufile, 0L, SEEK_END);
//
//  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"%s\""
//               " NumberOfComponents=\"3\" format=\"%s\">\n",
//           MANGLL_VTK_FLOAT_NAME, vectorName, MANGLL_VTK_FORMAT_STRING);
//
//  MANGLL_ASSERT (dim == 3);
//#ifdef MANGLL_VTK_ASCII
//  for (il = 0; il < Ntotal; ++il) {
//#ifdef MANGLL_VTK_DOUBLES
//    fprintf (vtufile, "     %24.16e %24.16e %24.16e\n",
//             v1d[vskip * il], v2d[vskip * il], v3d[vskip * il]);
//#else
//    fprintf (vtufile, "          %16.8e %16.8e %16.8e\n",
//             v1d[vskip * il], v2d[vskip * il], v3d[vskip * il]);
//#endif
//  }
//#else
//  float_data = MANGLL_ALLOC (MANGLL_VTK_FLOAT_TYPE, dim * Ntotal);
//  for (il = 0; il < Ntotal; ++il) {
//    float_data[3 * il + 0] = (MANGLL_VTK_FLOAT_TYPE) v1d[vskip * il];
//    float_data[3 * il + 1] = (MANGLL_VTK_FLOAT_TYPE) v2d[vskip * il];
//    float_data[3 * il + 2] = (MANGLL_VTK_FLOAT_TYPE) v3d[vskip * il];
//  }
//
//  fprintf (vtufile, "          ");
//  retval = mangll_vtk_write_binary (vtufile, (char *) float_data,
//                                    sizeof (*float_data) * dim * Ntotal);
//  fprintf (vtufile, "\n");
//  if (retval) {
//    MANGLL_LERROR ("mangll_vtk: Error encoding points\n");
//    fclose (vtufile);
//    return -1;
//  }
//  MANGLL_FREE (float_data);
//#endif
//  fprintf (vtufile, "        </DataArray>\n");
//
//  if (ferror (vtufile)) {
//    MANGLL_LERROR ("mangll_vtk: Error writing point vector\n");
//    fclose (vtufile);
//    return -1;
//  }
//  if (fclose (vtufile)) {
//    MANGLL_LERROR ("mangll_vtk: Error closing point vector\n");
//    return -1;
//  }
//  vtufile = NULL;
//
//  /* Only have the root write to the parallel vtk file */
//  if (procRank == 0) {
//    char                pvtufilename[BUFSIZ];
//    FILE               *pvtufile;
//    snprintf (pvtufilename, BUFSIZ, "%s.pvtu", baseName);
//
//    pvtufile = fopen (pvtufilename, "a");
//    if (!pvtufile) {
//      MANGLL_LERRORF ("Could not open %s for output!\n", vtufilename);
//      return -1;
//    }
//
//    fprintf (pvtufile, "      <PDataArray type=\"%s\" Name=\"%s\""
//                 " NumberOfComponents=\"3\" format=\"%s\"/>\n",
//             MANGLL_VTK_FLOAT_NAME, vectorName, MANGLL_VTK_FORMAT_STRING);
//
//    if (ferror (pvtufile)) {
//      MANGLL_LERROR ("mangll_vtk: Error writing parallel point vector\n");
//      fclose (pvtufile);
//      return -1;
//    }
//    if (fclose (pvtufile)) {
//      MANGLL_LERROR ("mangll_vtk: Error closing parallel point vector\n");
//      return -1;
//    }
//  }
//
//  return 0;
//}
//
//int
//mangll_vtk_write_footer (mangll_t * ma, const char *baseName)
//{
//  char                vtufilename[BUFSIZ];
//  int                 p;
//  int                 procRank = ma->mpirank;
//  int                 numProcs = ma->mpisize;
//  FILE               *vtufile;
//
//  /* Have each proc write to its own file */
//  snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", baseName, procRank);
//  vtufile = fopen (vtufilename, "a");
//  if (vtufile == NULL) {
//    MANGLL_LERRORF ("Could not open %s for output!\n", vtufilename);
//    return -1;
//  }
//
//  fprintf (vtufile, "      </PointData>\n");
//  fprintf (vtufile, "    </Piece>\n");
//  fprintf (vtufile, "  </UnstructuredGrid>\n");
//  fprintf (vtufile, "</VTKFile>\n");
//
//  if (ferror (vtufile)) {
//    MANGLL_LERROR ("mangll_vtk: Error writing footer\n");
//    fclose (vtufile);
//    return -1;
//  }
//  if (fclose (vtufile)) {
//    MANGLL_LERROR ("mangll_vtk: Error closing footer\n");
//    return -1;
//  }
//  vtufile = NULL;
//
//  /* Only have the root write to the parallel vtk file */
//  if (procRank == 0) {
//    char                pvtufilename[BUFSIZ];
//    FILE               *pvtufile;
//    snprintf (pvtufilename, BUFSIZ, "%s.pvtu", baseName);
//
//    pvtufile = fopen (pvtufilename, "a");
//    if (!pvtufile) {
//      MANGLL_LERRORF ("Could not open %s for output!\n", vtufilename);
//      return -1;
//    }
//
//    fprintf (pvtufile, "    </PPointData>\n");
//    for (p = 0; p < numProcs; ++p) {
//      fprintf (pvtufile, "    <Piece Source=\"%s_%04d.vtu\"/>\n", baseName,
//               p);
//    }
//    fprintf (pvtufile, "  </PUnstructuredGrid>\n");
//    fprintf (pvtufile, "</VTKFile>\n");
//
//    if (ferror (pvtufile)) {
//      MANGLL_LERROR ("mangll_vtk: Error writing parallel footer\n");
//      fclose (pvtufile);
//      return -1;
//    }
//    if (fclose (pvtufile)) {
//      MANGLL_LERROR ("mangll_vtk: Error closing parallel footer\n");
//      return -1;
//    }
//  }
//
//  return 0;
//}
//
