//
// Created by hari on 10/18/15.
//
// Mostly a port of VTK export functions in dendro
//

#ifndef DENDRO_VTK_H
#define DENDRO_VTK_H

#include <vector>
#include <TreeNode.h>
#include <string>


// base64 code from libb64
typedef enum
{
    step_A, step_B, step_C
} base64_encodestep;

typedef struct
{
    base64_encodestep step;
    char result;
    int stepcount;
} base64_encodestate;

void base64_init_encodestate(base64_encodestate* state_in);
char base64_encode_value(char value_in);
int base64_encode_block(const char* plaintext_in, int length_in, char* code_out, base64_encodestate* state_in);
int base64_encode_blockend(char* code_out, base64_encodestate* state_in);


/** This function writes numeric binary data in VTK base64 encoding.
 * \param vtkfile        Stream openened for writing.
 * \param numeric_data   A pointer to a numeric data array.
 * \param byte_length    The length of the data array in bytes.
 * \return               Returns 0 on success, -1 on file error.
 */
int                 dendro_vtk_write_binary (FILE * vtkfile, char *numeric_data,
                                         size_t byte_length);

/** This will write out the MPI rank in VTK format.
 *
 * This is a convenience function for the special
 * case of writing out the MPI rank only.  Note this
 * function will abort if there is a file error.
 *
 *  \param tree       The octree to be output.
 *  \param baseName The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be baseName_procNum.vtu).
 */
void                dendro_vtk_write_file (std::vector<ot::TreeNode> &tree,
                                           const char *baseName);

/** This will write the header of the vtu file.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of
 * fields.  The calling sequence would be something like
 *
 * \begincode
 * dendro_vtk_write_header(tree, "output", "Field1,Field2", "V", 1);
 * dendro_vtk_write_point_scalar(tree, "output", field1, 1, "Field1");
 * dendro_vtk_write_point_scalar(tree, "output", field2, 1, "Field2");
 * dendro_vtk_write_point_vector(tree, "output", v, 3, "V");
 * dendro_vtk_write_footer(tree, "output");
 * \endcode
 *
 *  \param tree         The tree to be outputted.
 *  \param baseName     The first part of the name which will have
 *                      the proc number appended to it (i.e., the
 *                      output file will be baseName_procNum.vtu).
 *  \param pointscalars A comma separated string of the scalar
 *                      fields that will be outputed at the nodes
 *                      of the mesh.  This can be set to \c NULL
 *                      if there are not any scalars.
 *  \param pointvectors A comma separated string of the vector
 *                      fields that will be outputed at the nodes
 *                      of the mesh.  This can be set to \c NULL
 *                      if there are not any scalars.
 *  \param write_rank   Boolean to determine if the MPI rank should be output.
 **
 *  \return             This returns 0 if no error and -1 if there is an error.
 */
int                 dendro_vtk_write_header (std::vector<ot::TreeNode> &tree,
                                             const char *baseName,
                                             const char *pointscalars,
                                             const char *pointvectors,
                                             int write_rank);

/** Extend dendro_vtk_write_header by writing wrapped mpi ranks.
 * \param [in] rank_wrap        If > 0, will wrap mpi ranks by modulo.
 */
int                 dendro_vtk_write_header_wrap (std::vector<ot::TreeNode>& tree,
                                                  const char *baseName,
                                                  const char *pointscalars,
                                                  const char *pointvectors,
                                                  int write_rank,
                                                  int rank_wrap);

/** Extend dendro_vtk_write_header_wrap by optionally writing EtoV.
 * \param [in] write_tag        If true, write the EtoV data to cells.
 */
int                 dendro_vtk_write_header_all (std::vector<ot::TreeNode>& tree,
                                                 const char *baseName,
                                                 const char *pointscalars,
                                                 const char *pointvectors,
                                                 int write_rank,
                                                 int rank_wrap,
                                                 int write_tag);

/** This will write a scalar field to the vtu file.
 *
 * It is good practice to make sure that the scalar field also
 * exists in the comma separated string \a pointscalars passed
 * to \c dendro_vtk_write_header.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of
 * fields.  The calling sequence would be something like
 *
 * \begincode
 * dendro_vtk_write_header(tree,  "output", "Field1,Field2", "V", 1);
 * dendro_vtk_write_point_scalar(tree, "output", field1, 1, "Field1");
 * dendro_vtk_write_point_scalar(tree, "output", field2, 1, "Field2");
 * dendro_vtk_write_point_vector(tree, "output", v, 3, "V");
 * dendro_vtk_write_footer(tree, "output");
 * \endcode
 *
 *  \param tree       The octree to be outputted.
 *  \param baseName   The first part of the name which will have
 *                    the proc number appended to it (i.e., the
 *                    output file will be baseName_procNum.vtu).
 *  \param scalard    The finite element field that will be
 *                    written out.
 *  \param sskip      The number of entries to advance in the array
 *                    to get the next scalar field entry.
 *  \param scalarName The name of the scalar field.
 **
 *  \return           This returns 0 if no error and -1 if there is an error.
 */
int                 dendro_vtk_write_point_scalar (std::vector<ot::TreeNode> &tree,
                                                   const char *baseName,
                                                   double *scalard,
                                                   int sskip,
                                                   const char *scalarName);

/** This will write a vector field to the vtu file.
 *
 * It is good practice to make sure that the vector field also
 * exists in the comma separated string \a pointvectors passed
 * to \c dendro_vtk_write_header.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of
 * fields.  The calling sequence would be something like
 *
 * \begincode
 * dendro_vtk_write_header(tree,  "output", "Field1,Field2", "V", 1);
 * dendro_vtk_write_point_scalar(tree, "output", field1, 1, "Field1");
 * dendro_vtk_write_point_scalar(tree, "output", field2, 1, "Field2");
 * dendro_vtk_write_point_vector(tree, "output", v, 3, "V");
 * dendro_vtk_write_footer(tree, "output");
 * \endcode
 *
 *  \param ma         The octree to be outputted.
 *  \param baseName   The first part of the name which will have
 *                    the proc number appended to it (i.e., the
 *                    output file will be baseName_procNum.vtu).
 *  \param vd         The x- y- and z-components, of the finite element vector
 *                    that will be written out.
 *  \param vskip      The number of entries to advance in the array
 *                    to get the next set of vector entries.
 *  \param vectorName The name of the scalar field.
 *
 *  \return           This returns 0 if no error and -1 if there is an error.
 */
int                 dendro_vtk_write_point_vector (std::vector<ot::TreeNode> &tree,
                                                   const char *baseName,
                                                   double *vd,
                                                   int vskip,
                                                   const char *vectorName);

/** This will write a vector field to the vtu file.
 *
 * It is good practice to make sure that the vector field also
 * exists in the comma separated string \a pointvectors passed
 * to \c dendro_vtk_write_header.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of
 * fields.  The calling sequence would be something like
 *
 * \begincode
 * dendro_vtk_write_header(tree,  "output", "Field1,Field2", "V", 1);
 * dendro_vtk_write_point_scalar(tree, "output", field1, 1, "Field1");
 * dendro_vtk_write_point_scalar(tree, "output", field2, 1, "Field2");
 * dendro_vtk_write_point_vector_123(tree, "output", v1, v2, v3, "V");
 * dendro_vtk_write_footer(tree, "output");
 * \endcode
 *
 *  \param ma         The octree to be outputted.
 *  \param baseName   The first part of the name which will have
 *                    the proc number appended to it (i.e., the
 *                    output file will be baseName_procNum.vtu).
 *  \param v1d        The x-components, of the finite element vector
 *                    that will be written out.
 *  \param v2d        The y-components, of the finite element vector
 *                    that will be written out.
 *  \param v3d        The z-components, of the finite element vector
 *                    that will be written out.
 *  \param vskip      The number of entries to advance in the array
 *                    to get the next set of vector entries.
 *  \param vectorName The name of the scalar field.
 *
 *  \return           This returns 0 if no error and -1 if there is an error.
 */
int                 dendro_vtk_write_point_vector_123 (std::vector<ot::TreeNode> &tree,
                                                       const char *baseName,
                                                       double *v1d,
                                                       double *v2d,
                                                       double *v3d,
                                                       int vskip,
                                                       const char
                                                       *vectorName);

/** This will write the footer of the vtu file.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of
 * fields.  To write out two fields the
 * calling sequence would be something like
 *
 * \begincode
 * dendro_vtk_write_header(tree, "output", "Field1,Field2", "V", 1);
 * dendro_vtk_write_point_scalar(tree, "output", field1, 1, "Field1");
 * dendro_vtk_write_point_scalar(tree, "output", field2, 1, "Field2");
 * dendro_vtk_write_point_vector(tree, "output", v, 3, "V");
 * dendro_vtk_write_footer(tree, "output");
 * \endcode
 *
 *  \param ma       The dendro to be outputted.
 *  \param baseName The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be baseName_procNum.vtu).
 *
 *  \return         This returns 0 if no error and -1 if there is an error.
 */
int                 dendro_vtk_write_footer (std::vector<ot::TreeNode> &tree, const char *baseName);

#endif //DENDRO_DENDRO_VTK_H
