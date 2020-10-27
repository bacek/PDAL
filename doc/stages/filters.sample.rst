.. _filters.sample:

filters.sample
===============================================================================

The **Sample Filter** performs Poisson sampling of the input ``PointView``. The
The practice of performing Poisson sampling via "Dart Throwing" was introduced
in the mid-1980's by [Cook1986]_ and [Dippe1985]_, and has been applied to
point clouds in other software [Mesh2009]_.

Our implementation of Poisson sampling is made streamable by voxelizing the
space and only adding points to the output ``PointView`` if they do not violate
the minimum distance criterion (as specified by ``radius``). The voxelization
allows several optimizations, first by checking for existing points within the
same voxel as the point under consideration, which are mostly likely to
violate the minimum distance criterion. Furthermore, we can easily visit
neighboring voxels (limiting the search to those that are populated) without
the need to create a KD-tree from the entire input ``PointView`` first and
performing costly spatial searches.

.. seealso::

    :ref:`filters.decimation`, :ref:`filters.fps`,
    :ref:`filters.relaxationdartthrowing`,
    :ref:`filters.voxelcenternearestneighbor`,
    :ref:`filters.voxelcentroidnearestneighbor`, and :ref:`voxeldownsize` also
    perform decimation.

.. note::

    Starting with PDAL v2.3, the ``filters.sample`` now supports streaming
    mode. As a result, there is no longer an option to ``shuffle`` points (or
    to provide a ``seed`` for the shuffle).

.. note::

    Starting with PDAL v2.3, a ``cell`` option has been added that works with
    the existing ``radius``. The user must provide one or the other, but not
    both. The provided option will be used to automatically compute the other.
    The relationship between ``cell`` and ``radius`` is such that the
    ``radius`` defines the radius of a sphere that circumscribes a voxel with
    edge length defined by ``cell``.

.. embed::

.. streamable::

Options
-------------------------------------------------------------------------------

cell
  Voxel cell size. If ``radius`` is set, ``cell`` is automatically computed
  such that the cell is circumscribed by the sphere defined by ``radius``.

radius
  Minimum distance between samples. If ``cell`` is set, ``radius`` is
  automatically computed to defined a sphere that circumscribes the voxel cell.
  Whether specified or derived, ``radius`` defines the minimum allowable
  distance between points.

.. include:: filter_opts.rst

