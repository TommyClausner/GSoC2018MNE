import numpy as np
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.metrics import CCMetric
import os.path
import imageio
import os
from PIL import Image

from dipy.viz import regtools


def exampleDipy():
    import ssl
    if hasattr(ssl, '_create_unverified_context'):
        ssl._create_default_https_context = ssl._create_unverified_context
    from dipy.data import fetch_stanford_hardi, read_stanford_hardi
    fetch_stanford_hardi()
    nib_stanford, gtab_stanford = read_stanford_hardi()
    stanford_b0 = np.squeeze(nib_stanford.get_data())[..., 0]

    from dipy.data.fetcher import fetch_syn_data, read_syn_data
    fetch_syn_data()
    nib_syn_t1, nib_syn_b0 = read_syn_data()
    syn_b0 = np.array(nib_syn_b0.get_data())

    from dipy.segment.mask import median_otsu

    stanford_b0_masked, stanford_b0_mask = median_otsu(stanford_b0, 4, 4)
    syn_b0_masked, syn_b0_mask = median_otsu(syn_b0, 4, 4)

    static = stanford_b0_masked
    static_affine = nib_stanford.affine
    moving = syn_b0_masked
    moving_affine = nib_syn_b0.affine

    pre_align = np.array(
        [[1.02783543e+00, -4.83019053e-02, -6.07735639e-02, -2.57654118e+00],
         [4.34051706e-03, 9.41918267e-01, -2.66525861e-01, 3.23579799e+01],
         [5.34288908e-02, 2.90262026e-01, 9.80820307e-01, -1.46216651e+01],
         [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00]])

    from dipy.align.imaffine import AffineMap
    affine_map = AffineMap(pre_align,
                           static.shape, static_affine,
                           moving.shape, moving_affine)

    resampled = affine_map.transform(moving)

    metric = CCMetric(3)

    level_iters = [10, 10, 5]
    sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)

    mapping = sdr.optimize(static, moving, static_affine, moving_affine,
                           pre_align)

    warped_moving = mapping.transform(moving)

    for slice in range(41 - 12, 41 + 13):
        regtools.overlay_slices(static, resampled, slice, 1, 'Static',
                                'Pre Moving',
                                'GIFexample1/' + str(slice) + 'T1pre.png')
        regtools.overlay_slices(static, warped_moving, slice, 1, 'Static',
                                'Post moving',
                                'GIFexample1/' + str(slice) + 'T1post.png')


def abstractExample():
    numBlobs = 4

    shape = np.array([81, 106, 76])
    static = np.zeros(shape)
    resampled = np.zeros(shape)

    size_blob = 3

    for blob in range(numBlobs):
        jitter = np.random.randint(5)

        randIndX = 41 + 2 * blob
        randIndY = np.random.randint(shape[1] - 2 * size_blob)
        randIndZ = np.random.randint(shape[2] - 2 * size_blob)

        static[randIndX - size_blob:randIndX + size_blob,
        randIndY - size_blob:randIndY + size_blob,
        randIndZ - size_blob:randIndZ + size_blob] = 1

        randIndY = randIndY + jitter
        randIndZ = randIndZ + jitter + np.random.randint(-1, 2)

        resampled[randIndX - size_blob:randIndX + size_blob,
        randIndY - size_blob:randIndY + size_blob,
        randIndZ - size_blob:randIndZ + size_blob] = 1

    metric = CCMetric(3)

    level_iters = [10, 10, 5]
    sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)

    mapping = sdr.optimize(static, resampled)

    warped_moving = mapping.transform(resampled)

    for slice in range(41 - 3 * size_blob, 41 + 3 * size_blob + 1):
        regtools.overlay_slices(static, resampled, slice, 0, 'Static',
                                'Pre Moving',
                                'GIFexample1/' + str(slice) + 'preMov.png')

        regtools.overlay_slices(static, warped_moving, slice, 0, 'Static',
                                'Post moving',
                                'GIFexample1/' + str(slice) + 'postMov.png')


def makeGIF(searchString):
    images = []
    for filename in sorted(os.listdir(os.getcwd() + '/GIFexample1')):
        if 'png' in filename and searchString in filename:
            images.append(
                imageio.imread(os.getcwd() + '/GIFexample1/' + filename))

    imageio.mimsave(os.getcwd() + '/GIFexample1/' + searchString + '.gif',
                    images)


def combineImages(searchStringA, searchStringB, newName):
    for filename in sorted(os.listdir(os.getcwd() + '/GIFexample1')):
        if 'png' in filename and searchStringA in filename:
            filename2 = filename.replace(searchStringA, searchStringB)
            images = map(Image.open, [os.getcwd() + '/GIFexample1/' + filename,
                                      os.getcwd() + '/GIFexample1/' + filename2])
            widths, heights = zip(*(i.size for i in images))

            total_width = max(widths)
            max_height = sum(heights)

            new_im = Image.new('RGB', (total_width, max_height))

            y_offset = 0
            for im in images:
                new_im.paste(im, (0, y_offset))
                y_offset += im.size[1]

            new_im.save(os.getcwd() + '/GIFexample1/' + filename[
                                                        0:2] + newName + '.png')


exampleDipy()
abstractExample()
combineImages('preMov', 'postMov', 'prepostMov')
combineImages('T1pre', 'T1post', 'T1prepost')
makeGIF('prepostMov')
makeGIF('T1prepost')
