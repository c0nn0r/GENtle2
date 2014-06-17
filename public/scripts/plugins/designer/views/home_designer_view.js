/**
@module Designer
@submodule Views
@class HomeDesignerView
**/
define(function(require) {
  var Backbone    = require('backbone.mixed'),
      template    = require('hbars!../templates/home_designer_view'),
      Sequence    = require('sequence/models/sequence'),
      Gentle      = require('gentle')(),
      HomeDesignerView;

  HomeDesignerView = Backbone.View.extend({
    manage: true,
    template: template,
    className: 'home-designer',

    events: {
      'submit .home-designer-form': 'createNewSequence',
    },

    createNewSequence: function(event) {
      event.preventDefault();
      var $form     = this.$('form').first(),
          name      = $form.find('[name=name]').val() || 'Unnamed',
          sequence  = new Sequence({
            name: name,
            sequence: '',
            displaySettings: {
              primaryView: 'designer'
            }
          });

      Gentle.addSequencesAndNavigate([sequence]);
    }
  });

  return HomeDesignerView;
});