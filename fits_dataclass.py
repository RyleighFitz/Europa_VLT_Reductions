from collections import defaultdict
import fitsio
import numpy as np
from typing import Dict, Any, List, Union, Callable, Set, Iterable, Iterator, Optional, Tuple
from dataclasses import dataclass


@dataclass
class Image:
    """ Representation of an image, along with metadata. """

    header: Dict[str, Any]
    _data: Optional[np.ndarray] = None
    _filename: Optional[str] = None

    @property
    def shape(self) -> Tuple[int]:
        return self.data.shape

    @property
    def data(self) -> np.ndarray:
        if self._data is not None:
            return self._data
        try:
            self._data = fitsio.read(self._filename,header=False)
            return self._data
        except Exception as e:
            raise ValueError(f"Unable to load FITS data: {e}")

    def metadata(self, *keys: str) -> Tuple[Any]:
        return tuple(self.header.get(key, None) for key in keys)

    @classmethod
    def load(cls, filename, lazy=True):
        if lazy:
            h = fitsio.read_header(filename)
            return Image(header=h, _data=None, _filename=filename)
        else:
            d,h = fitsio.read(filename, header=True)
            return Image(header=h, _data=d, _filename=filename)


class ImageSet:
    """ A queryable set of images. """

    def __init__(self, *imagesets: Union[Image, Iterable[Image]]):
        self.images = []
        for ims in imagesets:
            if isinstance(ims, Image):
                self.images.append(ims)
            else:
                for im in ims:
                    self.images.append(im)

    def __iter__(self) -> Iterator[Image]:
        return iter(self.images)

    def __getitem__(self, item: int) -> Image:
        return self.images[item]

    def __len__(self) -> int:
        return len(self.images)

    def query(self, **kwargs) -> 'ImageSet':
        """ Get a subset of images with metadata satisfying certain conditions. """
        def matches(image: Image) -> bool:
            return image.metadata(*kwargs.keys()) == tuple(kwargs.values())
        return ImageSet([img for img in self if matches(img)])

    def __repr__(self):
        return f"<ImageSet of size {len(self)}>"


@dataclass
class Transform:
    """ A transformation mapping an ImageSet to an ImageSet. """
    transform_op: Callable[..., Union[Image, Iterable[Image]]]
    partition_keys: Optional[List[str]] = None

    def __call__(self, images: ImageSet) -> ImageSet:
        def partition_value(image: Image):
            if self.partition_keys is not None:
                return image.metadata(*self.partition_keys)
            else:
                return id(image)

        partition = defaultdict(list)
        for image in images:
            partition[partition_value(image)].append(image)
            
        transformed = [img for _images in partition.values()
                           for img in list(self.transform_op(*_images))]
        
        result = ImageSet(transformed)

        return result