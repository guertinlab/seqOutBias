/// Profile proc macro

extern crate proc_macro;
extern crate proc_macro2;
extern crate quote;
extern crate darling;
extern crate syn;

use proc_macro2::TokenStream;
use quote::{quote, ToTokens};
use syn::{Generics, Ident, parse_macro_input, DeriveInput};
use darling::{ast::Data, FromDeriveInput, FromField, FromMeta};


#[derive(FromDeriveInput)]
#[darling(supports(struct_named))]
pub(crate) struct Receiver {
  ident: Ident,
  generics: Generics,
  data: Data<(), ReceiverField>,
}

impl Receiver {
  fn fields_to_emit(&self) -> Vec<(String, syn::Type)> {
    self.data
      .as_ref()
      .take_struct()
      .expect("Profile only supports named structs")
      .into_iter()
      .filter(|field| !field.skip )
      .map(|field| (field.name(),field.ty.clone()))
      .filter( |pair| pair.0.starts_with( "flag_" ) && pair.0 != "flag_profile" )
      .collect()
  }
}

impl ToTokens for Receiver {
  fn to_tokens(&self, tokens: &mut TokenStream) {
    let ident = &self.ident;
    let (impl_generics, ty_generics, where_clause) = self.generics.split_for_impl();
    let fields = self.fields_to_emit();
    let names: Vec<String> = fields.iter().map( |field| field.0[5..].to_string() ).collect();
    let fields_ident: Vec<Ident> = fields.iter().map( |field| Ident::from_string(&field.0).unwrap() ).collect();
    let accessors: Vec<TokenStream> = fields
      .iter()
      .map( | field | {
        let ftype = &field.1;
        let stype = quote!{#ftype}.to_string();

        match stype.as_str() {
          "Option < String >" => quote!{ .as_str().map(|s| s.to_string() ) },
          "u8"|"u16"|"u32" => quote!{ .as_integer().unwrap() as #ftype },
          "u64" => quote!{ .as_integer().unwrap() },
          "bool" => quote!{ .as_bool().unwrap() },
          _ => { quote!{ } }
        }
      })
      .collect();

    tokens.extend(quote! {
            #[automatically_derived]
            impl #impl_generics #ident #ty_generics #where_clause {

                pub fn apply_profile( &mut self, value: &Value ) {
                  #( { let field = value.get(#names); if field.is_some() { println!( "# Loading flag '{}' from profile file.", #names ); self.#fields_ident = field.unwrap() #accessors; } } )*
                }
            }
        })
  }
}

#[derive(FromField)]
#[darling(attributes(field_names))]
struct ReceiverField {
  ident: Option<Ident>,
  ty:	syn::Type,
  #[darling(default)]
  skip: bool,
}

impl ReceiverField {
  fn name(&self) -> String {
    self.ident
      .as_ref()
      .expect("Profile only supports named fields")
      .to_string()
  }
}

#[proc_macro_derive(Profile)]
pub fn profile_derive(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
  Receiver::from_derive_input(&parse_macro_input!(input as DeriveInput))
    .map(|receiver| quote!(#receiver))
    .unwrap_or_else(|err| err.write_errors())
    .into()
}
