# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

# General purpose libraries
import numpy as np
from tqdm import tqdm
from pathlib import Path

# Read annotations
import json

# Create labels and masks
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from PIL import Image, ImageDraw
from skimage import draw as skimage_draw
from skimage import morphology
from scipy import ndimage
from utils_plot import cmap_create_random
from skimage.measure import regionprops
from skimage.io import imread, imsave
from skimage.segmentation import clear_border
from skimage import measure
import warnings

# ---------------------------------------------------------------------------
#  Label images
# ---------------------------------------------------------------------------


#  Function to create label images
def create_label_images(path_sample, sample_name, path_save, label_name,
                        annotation_file='annotation.json'):
    """ Function to recursively scan a folder for annotations. Will then
    create a label image from these annotations and save them to a new
    folder.

    Parameters
    ----------
    path_sample : pathlib Path object
        Path to be recursively scanned.
    sample_name : str
        Name of the sample. This is usually the name of the subfolder
        containing the annotation.
    path_save : pathlib Path object
        Folder were label image will be saved.
    label_name : str
        suffix that will be added to the sample name to store label image.
    annotation_file : str, optional
        Name of annotation file, by default 'annotation.json'.
    """

    # Geojson importer
    annotationsGeoJson = GeojsonImporter()

    # Color maps
    cmap_reds = cm.Reds
    cmap_reds.set_under('k', alpha=0)
    cmap_random = cmap_create_random()  

    # Read json annotation file
    file_json = path_sample / annotation_file
    if not file_json.is_file():
        print(f'Annotation not found: {file_json}')
        return

    annot_dict, roi_size_all, image_size = annotationsGeoJson.load(file_json)

    # Create binary masks of nuclei
    binaryMasks = BinaryMaskGenerator(image_size=(image_size[1], 
                                                  image_size[0]),
                                                  erode_size=10,
                                                  obj_size_rem=500,
                                                  save_indiv=True)

    print('Creating binary masks ... this can take a little bit of time')
    mask_dict = binaryMasks.generate(annot_dict)

    # Remove overlap (can occur from manual corrections)

    mask_overlap = mask_dict['fill_indiv'].sum(axis=2)
    mask_overlap[mask_overlap < 2] = 0

    mask_labels = mask_dict['labels']
    mask_labels[mask_overlap > 0] = 0

    # Remove small objects
    mask_labels = morphology.remove_small_objects(mask_labels, min_size=100)
    props = regionprops(mask_labels)

    # Save label image
    name_save = path_save / f'{sample_name}__{label_name}_labels.png'
    imsave(name_save, mask_labels, check_contrast=False)

    # Plot summary plots for manual inspection
    fig, ax = plt.subplots(1, 2)
    fig.set_size_inches((10, 5))

    ax[0].imshow(mask_dict['edge'], cmap="Blues")
    ax[0].imshow(mask_dict['fill'], cmap="Blues", alpha=0.2)
    ax[0].imshow(mask_overlap, cmap=cmap_reds, clim=[0.99, 1.01])
    ax[0].get_xaxis().set_visible(False)
    ax[0].get_yaxis().set_visible(False)

    for prop in props:
        ax[0].text(prop.centroid[1],
                   prop.centroid[0],
                   f'{prop.label}',
                   fontsize=6, weight='bold',
                   verticalalignment='center', horizontalalignment='center')

    ax[1].imshow(mask_labels, cmap=cmap_random)
    ax[1].get_xaxis().set_visible(False)
    ax[1].get_yaxis().set_visible(False)

    ax[0].set_title(f'Segmentation results with overlap in red')
    ax[1].set_title(f'Label images with overlap removed')

    plt.tight_layout()
    name_save = path_save / f'{sample_name}__{label_name}_info.png'
    plt.savefig(name_save, dpi=300)

    plt.close()


def process_labels_closest(path_scan, suffix_label, suffix_save, truncate_distance=None):
    """
    Function to process label images. Will create two additional images
    1. Distance to the closted object, e.g. nuclei\
    2. Index of clostest object
    """

    # Recursively look for all segmentation masks
    for file_label in path_scan.rglob(f'*{suffix_label}'):
        print(f'Analyzing file {file_label}')

        # >>>> Nuclei label image
        img_labels = imread(file_label)
        props = regionprops(img_labels)
        nuc_labels = np.array([prop.label for prop in props])
        n_nuclei = len(nuc_labels)

        # Loop over all nuclei and create create distance map
        print(' Creating distance maps. This can take a while ...')
        dist_mat = np.zeros((img_labels.shape[0], img_labels.shape[1], n_nuclei))
        mask_fill_indiv = np.zeros((img_labels.shape[0], img_labels.shape[1], n_nuclei))

        for indx, obj_int in enumerate(tqdm(np.nditer(nuc_labels), total=n_nuclei)):

            # Create binary mask for current object and find contour
            img_label_loop = np.zeros((img_labels.shape[0], img_labels.shape[1]))
            img_label_loop[img_labels == obj_int] = 1
            mask_fill_indiv[:, :, indx] = img_label_loop

            dist_nuc = ndimage.distance_transform_edt(np.logical_not(img_label_loop))
            if truncate_distance:
                dist_nuc[dist_nuc > truncate_distance] = truncate_distance
            dist_mat[:, :, indx] = dist_nuc

        # >>> Condense distmap in two matrixes: index and distance to closest nucleus
        dist_nuc_ind_3D = np.argsort(dist_mat, axis=2)
        dist_nuc_dist_3D = np.take_along_axis(dist_mat, dist_nuc_ind_3D, axis=2)

        # For index: replace Python matrix index with actual index from label image
        name_save = Path(str(file_label).replace(suffix_label,suffix_save[0]))
        ind_nucleus_closest = np.zeros((img_labels.shape[0], img_labels.shape[1]))
        dist_nuc_ind_2D = np.copy(dist_nuc_ind_3D[:,:,0])

        for indx, obj_int in enumerate(np.nditer(nuc_labels)):
            ind_nucleus_closest[dist_nuc_ind_2D == indx] = obj_int

        if str(name_save) != str(file_label):
            imsave(name_save, ind_nucleus_closest.astype('uint16'), check_contrast=False)
        else:
            print(f'Name to save index matrix could not be established: {name_save}')

        # Save distances
        name_save = Path(str(file_label).replace(suffix_label,suffix_save[1]))
        if str(name_save) != str(file_label):
            imsave(name_save, dist_nuc_dist_3D[:, :, 0].astype('uint16'), check_contrast=False)
        else:
            print(f'Name to save index matrix could not be established: {name_save}')

# ---------------------------------------------------------------------------
# Classes to import annotations
# ---------------------------------------------------------------------------

class AnnotationImporter:
    """
    Base class to import manual annotations importer.
    """

    def load(self, path_open):
        """
        Load and annotations and return dictionary with annotations.
        """

        raise NotImplementedError('No load function defined for this class!')


class GeojsonImporter(AnnotationImporter):
    """
    Class to import manual annotations from GeoJson files.
    """

    def __init__(self, image_size=(2048, 2048)):
        """
        Initiate annotation dictionary.

        Args:
            image_size (tuple): size of image.

        """
        self.image_size = image_size

    def load(self, file_open):
        """
        Read folder content based on defined config.

        Args:
            file_open (string): file-name of annotation.

        Returns:
            annot_dict (dictionary): contains all annotated elements
            roi_size_all (list): contains size of each annotated element
        """

        with open(file_open, encoding='utf-8-sig') as fh:
            data_json = json.load(fh)

        # Overwrite default file size if bounding box is present
        if 'bbox' in data_json:
            self.image_size = (int(data_json['bbox'][2] - data_json['bbox'][0] + 1),
                               int(data_json['bbox'][3] - data_json['bbox'][1] + 1))

        # Loop over list and create simple dictionary & get size of annotations
        annot_dict = {}
        roi_size_all = {}

        skipped = []

        for feat_idx, feat in enumerate(data_json['features']):

            if feat['geometry']['type'] not in ['Polygon', 'LineString']:
                skipped.append(feat['geometry']['type'])
                continue

            key_annot = 'annot_' + str(feat_idx)
            annot_dict[key_annot] = {}
            annot_dict[key_annot]['type'] = feat['geometry']['type']
            annot_dict[key_annot]['pos'] = np.squeeze(np.asarray(feat['geometry']['coordinates']))
            annot_dict[key_annot]['properties'] = feat['properties']

            # Store size of regions
            if not (feat['properties']['label'] in roi_size_all):
                roi_size_all[feat['properties']['label']] = []

            roi_size_all[feat['properties']['label']].append(
                [annot_dict[key_annot]['pos'][:, 0].max() -
                 annot_dict[key_annot]['pos'][:, 0].min(),
                 annot_dict[key_annot]['pos'][:, 1].max()
                 - annot_dict[key_annot]['pos'][:, 1].min()])

        print('Skipped geometry type(s):', skipped)
        return annot_dict, roi_size_all, self.image_size


# ---------------------------------------------------------------------------
# Classes to generate masks
# ---------------------------------------------------------------------------


class MaskGenerator:
    """
    Base class for mask generators.
    """

    def __init__(self):
        pass

    def generate(self, annot_dict):
        """
        Generate the masks and return a dictionary.
        """
        raise NotImplementedError('No load function defined for this class!')

    def save(self, mask_dict, mask_key, file_name):
        """
        Save selected mask to a png file.

        Args:
            mask_dict (dictionary): dictionary with masks.
            mask_key (string): key for mask that should be saved.
            file_name (string): file-name for mask
        """

        if not (mask_key in mask_dict.keys()):
            print(f'Selected key ({mask_key})is not present in mask dictionary.')
            return

        # Save label - different labels are saved differently
        mask_save = mask_dict[mask_key]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            if mask_key is 'distance_map':
                imsave(file_name, mask_save)

            elif (mask_key is 'edge') or (mask_key is 'fill'):
                imsave(file_name, 255 * mask_save)

            elif mask_key is 'edge_weighted':
                mask_rescale = (mask_save - mask_save.min()) * 255 / (mask_save.max() - mask_save.min())
                mask_rescale = mask_rescale.astype('uint8')
                imsave(file_name, mask_rescale)

            else:
                imsave(file_name, mask_save.astype('float32'))


class BinaryMaskGenerator(MaskGenerator):
    """
    Create binary masks from dictionary with annotations. Depending on the
    annotation type, different masks are created. If masks are
        - polygons  : edge mask and a filled mask are created.
        - freelines  : only an edge mask is created.

    """

    def __init__(self, image_size=(2048, 2048), erode_size=5, obj_size_rem=500, save_indiv=False):
        self.erode_size = erode_size
        self.obj_size_rem = obj_size_rem
        self.save_indiv = save_indiv
        self.image_size = (image_size[1], image_size[0])

    def generate(self, annot_dict):
        """
        Create masks from annotation dictionary

        Args:
            annot_dict (dictionary): dictionary with annotations

        Returns:
            mask_dict (dictionary): dictionary with masks
        """

        # Get dimensions of image and created masks of same size
        # This we need to save somewhere (e.g. as part of the geojson file?)

        # Filled masks and edge mask for polygons
        mask_fill = np.zeros(self.image_size, dtype=np.uint8)
        mask_edge = np.zeros(self.image_size, dtype=np.uint8)
        mask_labels = np.zeros(self.image_size, dtype=np.uint16)

        rr_all = []
        cc_all = []

        if self.save_indiv is True:
            mask_edge_indiv = np.zeros(
                (self.image_size[0], self.image_size[1], len(annot_dict)), dtype=np.bool)
            mask_fill_indiv = np.zeros(
                (self.image_size[0], self.image_size[1], len(annot_dict)), dtype=np.bool)

        # Image used to draw lines - for edge mask for freelines
        im_freeline = Image.new('1', (self.image_size[1], self.image_size[0]), color=0)
        draw = ImageDraw.Draw(im_freeline)

        # Loop over all roi
        i_roi = 0
        #for roi_key, roi in annot_dict.items():
        for roi_key, roi in tqdm(annot_dict.items()):
            roi_pos = roi['pos']
            y_inverted = [self.image_size[0] - r - 1 for r in roi_pos[:, 1]]
            x = roi_pos[:, 0]

            # Check region type
            # freeline - line
            if roi['type'] in ['freeline', 'LineString']:

                # Loop over all pairs of points to draw the line

                for ind in range(roi_pos.shape[0] - 1):
                    line_pos = ((roi_pos[ind, 1], roi_pos[ind, 0], roi_pos[
                        ind + 1, 1], roi_pos[ind + 1, 0]))
                    draw.line(line_pos, fill=1, width=self.erode_size)

            # freehand - polygon
            elif roi['type'] in ['freehand', 'polygon', 'polyline', 'Polygon']:

                # Draw polygon
                rr, cc = skimage_draw.polygon(y_inverted, x)
                rr_perimeter, cc_perimeter = skimage_draw.polygon_perimeter(y_inverted, x, self.image_size)

                # Make sure it's not outside
                rr[rr < 0] = 0
                rr[rr > self.image_size[0]-1] = self.image_size[0]-1

                cc[cc < 0] = 0
                cc[cc > self.image_size[1]-1] = self.image_size[1]-1

                # Test if this region has already been added
                if any(np.array_equal(rr, rr_test) for rr_test in rr_all) and any(
                        np.array_equal(cc, cc_test) for cc_test in cc_all):
                    # print('Region #{} has already been used'.format(i + 1))
                    continue

                rr_all.append(rr)
                cc_all.append(cc)

                # Generate mask
                mask_fill_roi = np.zeros(self.image_size, dtype=np.uint8)
                mask_fill_roi[rr, cc] = 1
                mask_fill_roi[rr_perimeter, cc_perimeter] = 1

                # Erode to get cell edge - both arrays are boolean to be used as
                # index arrays later
                mask_fill_roi_erode = morphology.binary_erosion(
                    mask_fill_roi, np.ones((self.erode_size, self.erode_size)))
                mask_edge_roi = (mask_fill_roi.astype('int') -
                                 mask_fill_roi_erode.astype('int')).astype('bool')

                # Save array for mask and edge
                mask_fill[mask_fill_roi > 0] = 1
                mask_edge[mask_edge_roi] = 1
                mask_labels[mask_fill_roi > 0] = i_roi + 1

                if self.save_indiv is True:
                    mask_edge_indiv[:, :, i_roi] = mask_edge_roi.astype('bool')
                    mask_fill_indiv[:, :, i_roi] = mask_fill_roi.astype('bool')

                i_roi = i_roi + 1

            else:
                roi_type = roi['type']
                raise NotImplementedError(f'Mask for roi type "{roi_type}" can not be created')

        del draw

        # Convert mask from free-lines to numpy array
        mask_edge_freeline = np.asarray(im_freeline)
        mask_edge_freeline = mask_edge_freeline.astype('bool')

        # Post-processing of fill and edge mask - if defined
        mask_dict = {}
        if np.any(mask_fill):

            # (1) remove edges , (2) remove small  objects
            mask_fill = mask_fill & ~mask_edge
            mask_fill = morphology.remove_small_objects(
                mask_fill.astype('bool'), self.obj_size_rem)

            # For edge - consider also freeline edge mask
            mask_edge = mask_edge.astype('bool')
            mask_edge = np.logical_or(mask_edge, mask_edge_freeline)

            # Assign to dictionary for return
            mask_dict['edge'] = mask_edge
            mask_dict['fill'] = mask_fill.astype('bool')
            mask_dict['labels'] = mask_labels.astype('uint16')

            if self.save_indiv is True:
                mask_dict['edge_indiv'] = mask_edge_indiv
                mask_dict['fill_indiv'] = mask_fill_indiv
            else:
                mask_dict['edge_indiv'] = np.zeros(self.image_size + (1,), dtype=np.uint8)
                mask_dict['fill_indiv'] = np.zeros(self.image_size + (1,), dtype=np.uint8)

        # Only edge mask present
        elif np.any(mask_edge_freeline):
            mask_dict['edge'] = mask_edge_freeline
            mask_dict['fill'] = mask_fill.astype('bool')
            mask_dict['labels'] = mask_labels.astype('uint16')

            mask_dict['edge_indiv'] = np.zeros(self.image_size + (1,), dtype=np.uint8)
            mask_dict['fill_indiv'] = np.zeros(self.image_size + (1,), dtype=np.uint8)

        else:
            raise Exception('No mask has been created.')

        return mask_dict


class DistanceMapGenerator(MaskGenerator):
    """
    Create a distance transform from the edge. Stored as 16bit float, for
    display and saving converted to float32 (.astype('float32'))

    Requires that binary weights are calculated first, which is generated with
    the BinaryMaskGenerator (with the option flag_save_indiv=True).
    """

    def __init__(self, truncate_distance=None):
        self.truncate_distance = truncate_distance

    def generate(self, annot_dict, mask_dict):
        """
        Creates a distance map with truncated distance to the edge of the cell.

        Args:
            annot_dict (dictionary): dictionary with annotations
            mask_dict (dictionary): dictionary with masks containing at
                                    least binary masks

        Returns:
            mask_dict (dictionary): dictionary with additional weighted masks
        """

        mask_fill_indiv = mask_dict['fill_indiv']
        mask_edge_indiv = mask_dict['edge_indiv']
        dist_mat = np.ones(np.shape(mask_fill_indiv))

        for i_cell in range(mask_fill_indiv.shape[-1]):
            img_cell = mask_edge_indiv[
                       :, :, i_cell] + mask_fill_indiv[:, :, i_cell]

            dist_cell = ndimage.distance_transform_edt(img_cell)
            if self.truncate_distance:
                dist_cell[dist_cell >
                          self.truncate_distance] = self.truncate_distance
            dist_mat[:, :, i_cell] = dist_cell

        dist_map = np.sum(dist_mat, 2)

        # Note: saved as uint 16
        mask_dict['distance_map'] = dist_map.astype('uint16')

        return mask_dict



class BorderMaskGenerator(MaskGenerator):
    """
    https://github.com/selimsef/dsb2018_topcoders
    """

    def __init__(self, border_detection_threshold=6):
        self.border_detection_threshold = border_detection_threshold

    def generate(self, annot_dict, mask_dict):
        labels = mask_dict['labels']
        tmp = mask_dict['edge'] > 0
        tmp = morphology.dilation(tmp, morphology.square(self.border_detection_threshold))

        props = measure.regionprops(labels)
        msk0 = 255 * (labels > 0)
        msk0 = msk0.astype('uint8')

        msk1 = np.zeros_like(labels, dtype='bool')

        max_area = np.max([p.area for p in props])

        for y0 in range(labels.shape[0]):
            for x0 in range(labels.shape[1]):
                if not tmp[y0, x0]:
                    continue
                sz = self.border_detection_threshold

                uniq = np.unique(labels[max(0, y0 - sz):min(labels.shape[0], y0 + sz + 1),
                                 max(0, x0 - sz):min(labels.shape[1], x0 + sz + 1)])
                if len(uniq[uniq > 0]) > 1:
                    msk1[y0, x0] = True
                    msk0[y0, x0] = 0

        msk0 = 255 * (labels > 0)
        msk0 = msk0.astype('uint8')  # cell area
        msk1 = morphology.binary_closing(msk1)
        msk1 = 255 * msk1  # cell boundarys
        msk1 = msk1.astype('uint8')

        msk2 = np.zeros_like(labels, dtype='uint8')
        msk = np.stack((msk0, msk1, msk2))
        msk = np.rollaxis(msk, 0, 3)

        # Note: saved as float 16 - to plot has to be converted to float32
        # To be saved rescaled as 8 bit
        mask_dict['border_mask'] = msk.astype('float32')

        return mask_dict
