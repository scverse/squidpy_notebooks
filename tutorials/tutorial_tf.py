# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: sphinx
#       format_version: '1.1'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3.9.12 ('squidpy_3_9')
#     language: python
#     name: python3
# ---

"""
# Predict cluster labels spots using Tensorflow

In this tutorial, we show how you can use the `squidpy.im.ImageContainer` object to train a ResNet model to predict cluster labels of spots.

This is a general approach that can be easily extended to a variety of supervised, self-supervised or unsupervised tasks. We aim to highlight how the flexibility provided by the image container, and it's seamless integration with AnnData, makes it easy to interface your data with modern deep learning frameworks such as Tensorflow.

Furthermore, we show how you can leverage such a ResNet model to generate a new set of features that can provide useful insights on spots similarity based on image morphology. 

First, we'll load some libraries. Note that Tensorflow is not a dependency of Squidpy and you'd therefore have to install it separately in your conda environment. Have a look at [the Tensorflow installation instructions](https://www.tensorflow.org/install). This of course applies to any deep learning framework of your choice.
"""

import scanpy as sc
import squidpy as sq
from squidpy.im import ImageContainer
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import AnnData
from sklearn.model_selection import (
    train_test_split,
)  # we'll use this function to split our dataset in train and test set
import tensorflow as tf
from tensorflow.keras.layers.experimental import (
    preprocessing,
)  # let's use the new pre-processing layers for resizing and data augmentation tasks

sc.logging.print_header() # TODO: update Scanpy and Squidpy versions
print(f"squidpy=={sq.__version__}")
print(f"tensorflow=={tf.__version__}")

###############################################################################
# We will load the public data available in Squidpy.

adata = sq.datasets.visium_hne_adata()
img = sq.datasets.visium_hne_image()

###############################################################################
# ## Create train-test split
# We create a vector of our labels with which to train the classifier. In this case, we will train a classifier to predict cluster labels obtained from gene expression. We'll create a one-hot encoded array with the convenient function `tf.one_hot`. Furthermore, we'll split the vector indices to get a train and test set. Note that we specify the cluster labels as the `stratify` argument, to make sure that the cluster labels are balanced in each split. 

# get train,test split stratified by cluster labels
train_idx, test_idx = train_test_split(
    adata.obs_names.values,
    test_size=0.2,
    stratify=adata.obs["cluster"],
    shuffle=True,
    random_state=42,
)

""
print(
    f"Train set : \n {adata[train_idx, :].obs.cluster.value_counts()} \n \n Test set: \n {adata[test_idx, :].obs.cluster.value_counts()}"
)


###############################################################################
# ## Create datasets and train the model
# Next, we'll create a Tensorflow dataset which will be used as data loader for model training. A key aspect of this step is how the Image Container makes it easy to relate spots information to the underlying image.
# In particular, we will make use of `img.generate_spot_crops`, a method that creates a generator to crop the tissue image corresponding to each spot. 
# In just one line of code you can create this generator as well as specifying the size of the crops . You might want to increase the size to include some neighborhood morphology information. 
#
# We won't get too much in details of the additional arguments and steps related to the Tensorflow Dataset objects, you can familiarize yourself with Tensorflow datasets [here](https://www.tensorflow.org/api_docs/python/tf/data/Dataset).

def get_ohe(adata: AnnData, cluster_key: str, obs_names: np.ndarray):
    cluster_labels = adata[obs_names, :].obs["cluster"]
    classes = cluster_labels.unique().shape[0]
    cluster_map = {v: i for i, v in enumerate(cluster_labels.cat.categories.values)}
    labels = np.array([cluster_map[c] for c in cluster_labels], dtype=np.uint8)
    labels_ohe = tf.one_hot(labels, depth=classes, dtype=tf.float32)
    return labels_ohe


def create_dataset(
    adata: AnnData,
    img: ImageContainer,
    obs_names: np.ndarray,
    cluster_key: str,
    augment: bool,
    shuffle: bool,
):
    # image dataset
    spot_generator = img.generate_spot_crops(
        adata,
        obs_names=obs_names,  # this arguent specified the observations names
        scale=1.5,  # this argument specifies that we will consider some additional context under each spot. Scale=1 would crop the spot with exact coordinates
        as_array="image",  # this line specifies that we will crop from the "image" layer. You can specify multiple layers to obtain crops from multiple pre-processing steps.
        return_obs=False,
    )
    image_dataset = tf.data.Dataset.from_tensor_slices([x for x in spot_generator])

    # label dataset
    lab = get_ohe(adata, cluster_key, obs_names)
    lab_dataset = tf.data.Dataset.from_tensor_slices(lab)

    ds = tf.data.Dataset.zip((image_dataset, lab_dataset))

    if shuffle:  # if you want to shuffle the dataset during training
        ds = ds.shuffle(1000, reshuffle_each_iteration=True)
    ds = ds.batch(64)  # batch
    processing_layers = [
        preprocessing.Resizing(128, 128),
        preprocessing.Rescaling(1.0 / 255),
    ]
    augment_layers = [
        preprocessing.RandomFlip(),
        preprocessing.RandomContrast(0.8),
    ]
    if augment:  # if you want to augment the image crops during training
        processing_layers.extend(augment_layers)

    data_processing = tf.keras.Sequential(processing_layers)

    ds = ds.map(lambda x, y: (data_processing(x), y))  # add processing to dataset
    return ds


""
train_ds = create_dataset(adata, img, train_idx, "cluster", augment=True, shuffle=True)
test_ds = create_dataset(adata, img, test_idx, "cluster", augment=True, shuffle=True)

###############################################################################
# Here, we are actually instantiating the model. We'll use a pre-trained ResNet on ImageNet, and a dense layer for output. 

input_shape = (128, 128, 3)  # input shape
inputs = tf.keras.layers.Input(shape=input_shape)

# load Resnet with pre-trained imagenet weights
x = tf.keras.applications.ResNet50(
    weights="imagenet",
    include_top=False,
    input_shape=input_shape,
    classes=15,
    pooling="avg",
)(inputs)
outputs = tf.keras.layers.Dense(
    units=15,  # add output layer
)(x)
model = tf.keras.Model(inputs, outputs)  # create model
model.compile(
    optimizer=tf.keras.optimizers.Adam(learning_rate=1e-4),  # add optimizer
    loss=tf.keras.losses.CategoricalCrossentropy(from_logits=True),  # add loss
)

""
model.summary()

""
history = model.fit(
    train_ds,
    validation_data=test_ds,
    epochs=50,
    verbose=2,
)

###############################################################################
# We can plot training and test loss during training. Clearly it would benefit from some more fine-tuning :).

sns.lineplot(x=np.arange(50), y="loss", data=history.history)
sns.lineplot(x=np.arange(50), y="val_loss", data=history.history)

###############################################################################
# Calculate embedding and visualize results
# -----------------------------------------
#
# What we are actually interested in is the ResNet embedding values of the data after training. We expect that such an embedding contains relevant features of the image that can be used for downstream analysis such as clustering or integration with gene expression.
#
# For generating this embedding, we first create a new dataset, that contains the full list of spots, in the correct order and without augmentation.

full_ds = create_dataset(
    adata, img, adata.obs_names.values, "cluster", augment=False, shuffle=False
)

###############################################################################
# Then, we instantiate another model without the output layer, in order to get the final embedding layer.

model_embed = tf.keras.Model(inputs, x)
embedding = model_embed.predict(full_ds)

###############################################################################
# We can then save the embedding in a new AnnData, and copy over all the relevant metadata from the AnnData with gene expression counts...

adata_resnet = AnnData(embedding, obs=adata.obs.copy())
adata_resnet.obsm["spatial"] = adata.obsm["spatial"].copy()
adata_resnet.uns = adata.uns.copy()
adata_resnet

###############################################################################
# ... perform the standard clustering analysis.

sc.pp.scale(adata_resnet)
sc.pp.pca(adata_resnet)
sc.pp.neighbors(adata_resnet)
sc.tl.leiden(adata_resnet, key_added="resnet_embedding_cluster")
sc.tl.umap(adata_resnet)

###############################################################################
# Interestingly, it seems that despite the poor performance on the test set, the model has encoded some information relevant to separate spots from each other. The clustering annotation also resembles the original annotation based on gene expression similarity.

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.pl.umap(
    adata_resnet, color=["cluster", "resnet_embedding_cluster"], size=100, wspace=0.7
)

###############################################################################
# We can visualize the same information in spatial coordinates. Again some clusters seems to closely recapitulate the Hippocampus and Pyramidal layers clusters. It seems to have worked surprisingly well!

sq.pl.spatial_scatter(
    adata_resnet,
    color=["cluster", "resnet_embedding_cluster"],
    frameon=False,
    wspace=0.5,
)

###############################################################################
# An additional analysis could be to integrate information of both gene expression and the features learned by the ResNet classifier, in order to get a joint representation of both gene expression and image information. Such integration could be done for instance by concatenating the resulting PCA from the gene expression `adata` and the ResNet embedding `adata_resnet`. After concatenating the principal components, you could follow the usual steps of building a KNN graph and clustering with the leiden algorithm. 
#
# With this tutorial we have shown how to interface the Squidpy workflow with modern deep learning frameworks, and have inspired you with additional analysis that leverage several data modalities and powerful DL-based representations.
