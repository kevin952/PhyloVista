from django.urls import path
from .views import phylo_tree, home
urlpatterns = [
    path('', home, name='home'),
    path('phylo_tree/', phylo_tree, name='phylo_tree')
    # Define your app-specific URLs here
    # For example:
    # path('some-view/', views.some_view, name='some_view'),
]